# VIDRA: Variant-Informed Dose-Response Analysis

> Code for the manuscript: *Genetic dose-response modelling predicts drug mechanisms, dosing, and adverse events*

VIDRA is a statistical framework that models how genetic variation across the allele frequency and functional consequence spectrum informs gene-phenotype dose-response relationships. We integrated over 1.5 million associated variants from Open Targets data to derive more than 88,000 dose-response-like relationships, increasing gene-phenotype associations 3-fold when incorporating rare variants. We benchmarked VIDRA against known drug targets to derive a Therapeutic Potential Score, identifying over 2,000 high-potential targets, most of which are currently untargeted.

---

## Pipeline overview

The VIDRA v2 pipeline runs on **Google Cloud** in three sequential steps. Each step reads from and writes to `gs://vidra-2-0/`.

```
Step 1 — Data preparation      (PySpark on Dataproc Serverless)
   ↓
Step 2 — Variant annotation     (GCE VM: VEP 111 Docker + FoldX local lookup)
   ↓
Step 3 — Bayesian analysis      (PySpark + CmdStanPy on Dataproc Serverless)
```

### Step 1: Data preparation

**Script:** `scripts/pyspark_scripts/prepare_analysis_input.py`
**Runs on:** Dataproc Serverless (PySpark)
**Reads:** Open Targets data releases, AZ PheWAS CSV, gene mappings (all from `gs://vidra-2-0/raw_data/`)
**Writes:** `gs://vidra-2-0/vidra_analysis_ready/` (parquet partitioned by gene)

Ingests and harmonises four variant sources (AZ PheWAS rare variants, ClinVar rare variants, colocalised common variants, coding GWAS variants), joins burden tests, fills defaults, and writes analysis-ready parquet.

```bash
gcloud dataproc batches submit pyspark \
  scripts/pyspark_scripts/prepare_analysis_input.py \
  --project=open-targets-genetics-dev --region=europe-west1 \
  --deps-bucket=gs://vidra-2-0 \
  --container-image=europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \
  --properties=spark.sql.execution.arrow.pyspark.enabled=true \
  -- --bucket_name=vidra-2-0
```

### Step 2: Variant annotation

**Script:** `scripts/pyspark_scripts/annotate_variants_cli.py`
**Runs on:** GCE VM (`e2-standard-8`, 32 GB RAM, 250 GB SSD)
**Reads:** `gs://vidra-2-0/vidra_analysis_ready_manifest/` (unique variant IDs)
**Writes:** `gs://vidra-2-0/variant_annotations/<output_name>.parquet`

Annotates ~1.6M unique variants using:
- **VEP 111** (Docker: `ensemblorg/ensembl-vep:release_111.0`) with plugins: Blosum62, CADD, AlphaMissense, REVEL, PrimateAI, Conservation (GERP)
- **FoldX local lookup** for coding variants only: FoldX ΔΔG, pLDDT
- **VEP cache** for clinical significance

The VM is created and configured by `scripts/pyspark_scripts/run_annotation_vm.sh`:

```bash
# Create the VM and run annotation autonomously:
bash scripts/pyspark_scripts/run_annotation_vm.sh

# Or, with --no-delete to keep the VM for debugging:
bash scripts/pyspark_scripts/run_annotation_vm.sh --no-delete

# Isolated run with output suffix (reads/writes suffixed directories):
bash scripts/pyspark_scripts/run_annotation_vm.sh --suffix=_dev
```

To run the annotation script manually on the VM:

```bash
python3 /opt/vidra/annotate_variants_cli.py \
  --bucket_name vidra-2-0 \
  --vep_cache_dir /opt/vep/cache \
  --plugin_dir /opt/vep/plugins \
  --plugin_data_dir /opt/vep/plugin_data \
  --work_dir /tmp/vidra_annotation \
  --use_docker --docker_image ensemblorg/ensembl-vep:release_111.0 \
  --foldx_file /opt/vep/plugin_data/foldx_energy.csv.gz \
  --output_name variant_annotations.parquet \
  --threads 8 --buffer_size 5000

# To match a Step 1 suffix, add --output_suffix:
# ... --output_suffix _dev
```

#### Annotation output columns

| Column | Source | Transform |
|--------|--------|-----------|
| `as_blosum62` | VEP Blosum62 plugin | sigmoid 1/(1+exp(-x)), fillna(1) |
| `as_conservation` | VEP Conservation (GERP RS) | normalise [-12.36, 6.18] → [0,1], invert; NaN preserved for per-gene mean imputation in Step 3 |
| `as_sift` | VEP cache | raw score; NaN preserved for per-gene mean imputation in Step 3 |
| `as_polyphen` | VEP cache | raw score; NaN preserved; inverted (`1 - x`) in Step 3 |
| `as_cadd` | VEP CADD plugin | phred/50 clamped [0,1]; NaN preserved; inverted (`1 - x`) in Step 3 |
| `as_alphamissense` | VEP AlphaMissense plugin | raw am_pathogenicity; NaN preserved; inverted (`1 - x`) in Step 3 |
| `as_revel` | VEP REVEL plugin | raw score; NaN preserved; inverted (`1 - x`) in Step 3 |
| `as_primateai` | VEP PrimateAI plugin | raw score, fillna(0) |
| `as_loftool` | — | placeholder, always 0.0 |
| `as_plddt` | FoldX local lookup | plddt/100, fillna(0); inverted (`1 - x`) in Step 3 |
| `as_consequence` | VEP most_severe_consequence | ordinal 0–12 (VEP severity order); normalised (`1 - code/12`) in Step 3 |
| `as_clinicalSignificance` | VEP cache (`--check_existing`) | 17 ordered categories → integer codes 0–16 → MinMaxScale [0,1] |
| `foldxDdq_raw` | FoldX local lookup | raw ΔΔG, per-gene transform applied in Step 3 |

NaN preservation: missense-specific annotations (sift, polyphen, cadd, alphamissense, revel, conservation) output NaN for non-missense variants so that Step 3 can apply per-gene mean imputation.

### Step 3: Bayesian analysis

**Script:** `scripts/pyspark_scripts/run_bayesian_analysis.py`
**Runs on:** Dataproc Serverless (PySpark + CmdStanPy)
**Reads:** `gs://vidra-2-0/vidra_analysis_ready/` (Hive-partitioned by gene from Step 1), `gs://vidra-2-0/variant_annotations/` (38 MB lookup table from Step 2)
**Writes:** `gs://vidra-2-0/vidra_results/` (parquet partitioned by gene)

For each (gene, disease) pair:
1. Joins annotation lookup table (broadcast) onto analysis-ready data
2. Applies direction inversions (`1 - x` for revel, polyphen, cadd, alphamissense, plddt)
3. Normalises `as_consequence` (`1 - code/12`) and rescales `as_clinicalSignificance`
4. Computes per-gene FoldX transform (`as_foldx`)
5. Per-gene mean imputation for missense-specific annotations (revel, alphamissense, sift, polyphen, foldx, cadd, conservation) — NaN values are filled with the gene mean of observed values, then LoF variants are hardcoded to 0.0 (maximally damaging)
6. Deduplicates variants within (disease, source, QTL) groups
7. Runs single-variant or multi-variant Stan model via variational inference
8. Extracts posterior summaries (percentiles, posterior probabilities)

When using `--gene_list`, Spark uses an `.isin()` filter for guaranteed partition pruning,
so a 3-gene test reads only 3 partitions and requires minimal executors.

```bash 
# Production — all genes:
gcloud dataproc batches submit pyspark \
  scripts/pyspark_scripts/run_bayesian_analysis.py \
  --project=open-targets-genetics-dev --region=europe-west1 \
  --deps-bucket=gs://vidra-2-0 \
  --container-image=europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \
  --properties=spark.sql.execution.arrow.pyspark.enabled=true \
  -- --bucket_name=vidra-2-0

# Specific genes only — create a text file of ENSG IDs (one per line) and upload:
gsutil cp my_genes.txt gs://vidra-2-0/my_genes.txt
# Then run with --gene_list:
... -- --bucket_name=vidra-2-0 --gene_list=gs://vidra-2-0/my_genes.txt

# Specific genes + random padding (e.g. 3 target genes padded to 50 total):
... -- --bucket_name=vidra-2-0 --gene_list=gs://vidra-2-0/my_genes.txt --test_mode --test_genes 50

# Random subset only:
... -- --bucket_name=vidra-2-0 --test_mode --test_genes 50

# Isolated run with output suffix (must match Step 1 & 2 suffix):
... -- --bucket_name=vidra-2-0 --output_suffix=_dev
```

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `--bucket_name` | *(required)* | GCS bucket name (without `gs://`) |
| `--h1` | `0.1` | Hyperparameter h1 passed to Stan models |
| `--gene_list` | — | Path to a GCS text file of ENSG IDs (one per line) to restrict the analysis |
| `--test_mode` | off | Restrict to a random subset of genes. If `--gene_list` is also set, fills remaining slots randomly up to `--test_genes` |
| `--test_genes` | `50` | Total number of genes when `--test_mode` is active |
| `--output_suffix` | `""` | Suffix for input/output directory names (e.g. `_dev`). Must match the suffix used in Steps 1 & 2 |

---

## Repository structure

```
VIDRA/
├── README.md                          ← This file
├── stan_models/
│   ├── VIDRA.stan                     ← Multi-variant hierarchical model
│   └── VIDRA_single_variant.stan      ← Single-variant model
├── Docker_container_VIDRA_image/
│   ├── Dockerfile                     ← Dataproc Serverless image (Spark + CmdStanPy + Stan)
│   └── python_requirements.txt
├── scripts/
│   ├── pyspark_scripts/
│   │   ├── prepare_analysis_input.py  ← Step 1: data ingestion & preparation
│   │   ├── annotate_variants_cli.py   ← Step 2: VEP + FoldX local annotation
│   │   ├── run_bayesian_analysis.py   ← Step 3: Bayesian dose-response models
│   │   └── run_annotation_vm.sh       ← Creates & configures the GCE annotation VM
│   ├── setup/
│   │   ├── copy_opentargets_data.sh   ← One-off: copy OT data releases to GCS
│   │   └── download_cadd_to_gcs.sh    ← One-off: download CADD v1.7 (~80 GB) to GCS
│   ├── modules_py/                    ← Original VIDRA Python module (reference)
│   ├── ingest_data/                   ← Original per-source ingest scripts (reference)
│   ├── preStat_processing/            ← Original pre-processing scripts (reference)
│   ├── stat/                          ← Original Stan model estimation scripts (reference)
│   ├── post_stat/                     ← Result analysis and enrichment scripts
│   └── VIDRA_QC_analysis_scripts/     ← Power analysis scripts
└── Nextflow_pipeline/                 ← Original Nextflow pipeline (reference, replaced by PySpark)
```

### One-off setup scripts

These scripts only need to be run once to populate the GCS bucket with reference data:

```bash
# 1. Copy Open Targets data releases into gs://vidra-2-0/raw_data/open_targets/
bash scripts/setup/copy_opentargets_data.sh

# 2. Download CADD v1.7 data (~80 GB) to gs://vidra-2-0/raw_data/CADD/
bash scripts/setup/download_cadd_to_gcs.sh
```

### Python environments

**Steps 1 & 3** (Dataproc Serverless) use the pre-built Docker image with all dependencies installed: `europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1`

To build and push the image (requires Docker with `buildx` and Artifact Registry write access):

```bash
# Authenticate Docker to Artifact Registry first:
gcloud auth configure-docker europe-west1-docker.pkg.dev

# Build and push (must be linux/amd64 for Dataproc workers):
cd /path/to/VIDRA
docker buildx build \
  --platform linux/amd64 \
  -t europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \
  -f Docker_container_VIDRA_image/Dockerfile \
  --push \
  .
```

**Step 2** (GCE annotation VM) dynamically creates a Python venv during startup (see `run_annotation_vm.sh` Phase 4).

---

## Stan models

Both models are in `stan_models/` and are pre-compiled into the Dataproc Docker image.

- **`VIDRA.stan`** — Multi-variant hierarchical model. Estimates `slope_random[1..5]` per variant source and a `meta_slope` across sources. Currently uses 4 active annotation features as predictors: `as_revel`, `as_cadd`, `as_alphamissense` (protein_prior) and `as_clinicalSignificance` (disease_prior). Several additional annotations are declared but commented out in the Stan code (conservation, sift, polyphen, blosum62, foldx, consequence, plddt, primateai). Key parameters: `bO`/`bOse` (burden test log-OR and its SE); `yOR`/`yORse` (variant effect log-OR and SE); `xcest`/`yORest` have weak unconditional priors `normal(0, 10)` to regularise ClinVar-only groups where these are otherwise unconstrained.
- **`VIDRA_single_variant.stan`** — Single-variant model. Declares 7 annotation features in the data block (`as_conservation`, `as_sift`, `as_polyphen`, `as_cadd`, `as_alphamissense`, `as_clinicalSignificance`, `as_primateai`) but only 5 are actively used in the model: `as_conservation`, `as_cadd`, `as_alphamissense` (protein_prior) and `as_clinicalSignificance`, `as_primateai` (disease_prior). `as_sift` and `as_polyphen` are declared but unused.

---

## GCS bucket layout

```
gs://vidra-2-0/
├── raw_data/                          ← Input data (OT releases, AZ PheWAS, plugin data)
│   ├── open_targets/                  ← Coloc, GWAS, ClinVar, burden, variants, studies
│   ├── CADD/                          ← whole_genome_SNVs.tsv.gz (~80 GB)
│   ├── AlphaMissense_hg38.tsv.gz
│   ├── new_tabbed_revel.tsv.gz
│   ├── PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz
│   └── gerp_conservation_scores.homo_sapiens.GRCh38.bw
├── vidra_analysis_ready/              ← Step 1 output (1 parquet file per gene partition)
├── vidra_analysis_ready_manifest/     ← Unique variant IDs for annotation
├── variant_annotations/               ← Step 2 output (38 MB parquet lookup table)
├── vidra_results/                     ← Step 3 output (posterior summaries)
└── scripts/                           ← Deployed copies of pipeline scripts
```

---

## Infrastructure

| Component | Spec |
|-----------|------|
| **GCP project** | `open-targets-genetics-dev` |
| **Region** | `europe-west1` |
| **Dataproc Serverless** | Steps 1 & 3 (PySpark) |
| **Docker image** | `europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1` |
| **Annotation VM** | `vidra-annotation-vm`, `e2-standard-8` (8 vCPU, 32 GB), 250 GB SSD |
| **VEP Docker image** | `ensemblorg/ensembl-vep:release_111.0` |
| **GCS bucket** | `gs://vidra-2-0/` |

---

## Legacy code

The `Nextflow_pipeline/`, `scripts/ingest_data/`, `scripts/preStat_processing/`, `scripts/stat/`, and `scripts/modules_py/` directories contain the original pipeline code from the manuscript. These are retained for reference but are **not used** by the current PySpark pipeline.
