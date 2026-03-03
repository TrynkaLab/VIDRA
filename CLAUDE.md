# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

VIDRA (Variant-Informed Dose-Response Analysis) models how genetic variation informs gene-phenotype dose-response relationships. It runs on Google Cloud in three sequential steps, processing ~1.5M variants from Open Targets to produce 88,000+ gene-disease dose-response estimates.

## Pipeline architecture

```
Step 1 — scripts/pyspark_scripts/prepare_analysis_input.py   (Dataproc Serverless)
   ↓   reads: gs://vidra-2-0/raw_data/
   ↓   writes: gs://vidra-2-0/vidra_analysis_ready/ (parquet, partitioned by gene)
Step 2 — scripts/pyspark_scripts/annotate_variants_cli.py    (GCE VM: VEP 111 + FoldX)
   ↓   reads: gs://vidra-2-0/vidra_analysis_ready_manifest/
   ↓   writes: gs://vidra-2-0/variant_annotations/ (38 MB parquet lookup)
Step 3 — scripts/pyspark_scripts/run_bayesian_analysis.py    (Dataproc Serverless)
       reads: gs://vidra-2-0/vidra_analysis_ready/ + variant_annotations/
       writes: gs://vidra-2-0/vidra_results/
```

Steps 1 and 3 use a pre-built Docker image; Step 2 runs on a GCE VM created by `scripts/pyspark_scripts/run_annotation_vm.sh`.

## Running the pipeline

### Step 1: Data preparation (full run)
```bash
gcloud dataproc batches submit pyspark \
  scripts/pyspark_scripts/prepare_analysis_input.py \
  --project=open-targets-genetics-dev --region=europe-west1 \
  --deps-bucket=gs://vidra-2-0 \
  --container-image=europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \
  --properties=spark.sql.execution.arrow.pyspark.enabled=true \
  -- --bucket_name=vidra-2-0

# Test mode (200 random genes):
-- --bucket_name=vidra-2-0 --test_mode --test_genes 200
```

### Step 2: Variant annotation
```bash
# Creates VM, runs annotation, uploads to GCS, optionally self-deletes:
bash scripts/pyspark_scripts/run_annotation_vm.sh
bash scripts/pyspark_scripts/run_annotation_vm.sh --no-delete  # keep VM for debugging
```

### Step 3: Bayesian analysis
```bash
# Full run:
gcloud dataproc batches submit pyspark \
  scripts/pyspark_scripts/run_bayesian_analysis.py \
  --project=open-targets-genetics-dev --region=europe-west1 \
  --deps-bucket=gs://vidra-2-0 \
  --container-image=europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \
  --properties=spark.sql.execution.arrow.pyspark.enabled=true \
  -- --bucket_name=vidra-2-0

# Specific genes (upload a text file of ENSG IDs first):
gsutil cp vidra_test_genes.txt gs://vidra-2-0/my_genes.txt
-- --bucket_name=vidra-2-0 --gene_list=gs://vidra-2-0/my_genes.txt

# Specific genes + random padding to N total:
-- --bucket_name=vidra-2-0 --gene_list=gs://vidra-2-0/my_genes.txt --test_mode --test_genes 50
```

### Docker image (rebuild after changing Stan models or Dockerfile)
```bash
gcloud auth configure-docker europe-west1-docker.pkg.dev
docker buildx build \
  --platform linux/amd64 \
  -t europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \
  -f Docker_container_VIDRA_image/Dockerfile \
  --push .
```
**Important:** Stan models are compiled into the Docker image at build time (`/opt/vidra/stan_models/`). Any change to a `.stan` file requires rebuilding and pushing the image before re-running Step 3.

## Key data model

### Variant source encoding
- `GsourceLab=0`: Common QTL variants (colocalised eQTL/pQTL from OT)
- `GsourceLab=1`: AZ PheWAS rare variants
- `GsourceLab=2`: ClinVar rare variants
- `GsourceLab=3`: Coding GWAS variants (non-colocalising)
- `GqtlLab=0`: eQTL, `GqtlLab=1`: pQTL, `GqtlLab=2`: no QTL

### Stan slope_random mapping
`COMBINATION_SLOPE` in `run_bayesian_analysis.py` maps (GsourceLab, GqtlLab) pairs to `slope_random` indices:
- `slope_random[1]` = eQTL (0,0)
- `slope_random[2]` = pQTL (0,1)
- `slope_random[3]` = coding GWAS (3,2)
- `slope_random[4]` = AZ rare (1,2)
- `slope_random[5]` = ClinVar rare (2,2)

### Annotation convention
All protein function annotations follow: **0 = damaging/pathogenic, 1 = benign/tolerated**.
`INVERSION_COLS` in `run_bayesian_analysis.py` lists columns inverted (`1 - x`) in Step 3 because their raw scale is inverted (high = damaging): currently `as_revel`, `as_polyphen`, `as_cadd`, `as_alphamissense`, `as_plddt`.
`as_conservation` is inverted during Step 2 annotation (not Step 3).

### Default fill values
When effect-size data is missing, defaults match the original pipeline (`scripts/ingest_data/pre_processing_VIDRA_per_gene_pheno.py`):
- `xc=0.0, xcse=0.1` (no QTL effect)
- `yc=0.0, ycse=0.14` (no GWAS effect)
- `bO=0.0, bOse=2.0` (no burden signal)
- `as_clinicalSignificance=0.0` (not provided/benign) — the only annotation in Step 1 output

VEP-derived protein annotations (blosum62, sift, polyphen, cadd, revel, etc.) are NOT carried in Step 1 output. They come exclusively from the Step 2 annotation parquet, joined in Step 3.
In Step 2/3, most annotations default to 1.0 (benign). Exceptions: `as_primateai=0.0`, `as_clinicalSignificance=0.0`, `as_loftool=0.0`, `as_plddt=0.0` (in Step 2; 1.0 in Step 3 Spark fill).

### Per-gene mean imputation
Missense-specific annotations (`GENE_MEAN_IMPUTE_COLS`: `as_revel`, `as_alphamissense`, `as_sift`, `as_polyphen`, `as_foldx`, `as_cadd`, `as_conservation`) are left as NaN in Step 2 output for non-missense variants. In Step 3 `process_gene()`, NaN values are filled with the **per-gene mean** of observed values. If all variants in a gene have NaN for a column, the global `ANNOTATION_DEFAULTS` benign value is used as fallback. Loss-of-function variants (`stop_gained`, `start_lost`) are hardcoded to 0.0 (maximally damaging) after imputation. This matches the original pipeline (`pre_processing_VIDRA_per_gene_pheno.py` lines 440–449).

## Stan models

Both models are in `stan_models/`. Changes require rebuilding the Docker image.

- **`VIDRA.stan`** — Multi-variant hierarchical model used when a (gene, disease) pair has >1 unique variant. Estimates `slope_random[1..5]` per source and a hierarchical `slope`. Currently uses 4 active annotations: `as_revel`, `as_cadd`, `as_alphamissense` (protein_prior) and `as_clinicalSignificance` (disease_prior); others are commented out. `xcest` has a tight `normal(0, 0.2)` prior; `intercept` has a weak `normal(0, 10)` prior. ADVI: tries fullrank first, falls back to meanfield for large N (>~100 variants, where fullrank's O(P²) covariance over P=4N+8 parameters becomes ill-conditioned).
- **`VIDRA_single_variant.stan`** — Single-variant model. Declares 7 annotations in the data block but only 5 are actively used in the model: `as_conservation`, `as_cadd`, `as_alphamissense` (protein_prior) and `as_clinicalSignificance`, `as_primateai` (disease_prior). `as_sift` and `as_polyphen` are declared but unused. The `slope ~ normal(yOR/xc, ...)` division only executes inside `if (numG1 == 0)`, so `xc=0.0` for non-QTL variants is safe — Stan never evaluates that branch for them.

## Code layout

```
scripts/pyspark_scripts/   ← PRODUCTION pipeline (the three steps + VM script)
stan_models/               ← Stan model definitions
Docker_container_VIDRA_image/  ← Dockerfile + Python deps
scripts/post_stat/         ← Post-analysis enrichment and validation scripts
scripts/setup/             ← One-off GCS setup scripts (run once)
Nextflow_pipeline/         ← LEGACY reference (replaced by PySpark)
scripts/ingest_data/       ← LEGACY reference
scripts/stat/              ← LEGACY reference (original Stan estimation scripts)
scripts/modules_py/        ← LEGACY reference
scripts/preStat_processing/ ← LEGACY reference
```

Legacy directories are retained as reference for the original Nextflow/BigQuery pipeline described in the manuscript. Do not use them for the current pipeline.

## Infrastructure
- **GCP project**: `open-targets-genetics-dev`, **region**: `europe-west1`
- **GCS bucket**: `gs://vidra-2-0/`
- **Dataproc image**: `europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1`
- **Annotation VM**: `vidra-annotation-vm`, `e2-standard-8` (8 vCPU, 32 GB RAM, 250 GB SSD)
- **VEP**: `ensemblorg/ensembl-vep:release_111.0`
- **CmdStanPy**: 1.2.4, **CmdStan**: 2.34.1

## Known design decisions and gotchas

- **`bOse` is SE, not SD**: The burden test SE is `(log(UCI) - log(LCI)) / 3.92` (log-scale). The original pipeline mistakenly used raw OR scale and multiplied by `sqrt(N)`.
- **AZ SE fix**: Uses `(log(UCI) - log(LCI)) / 3.92` instead of the original `sqrt(50) * (UCI - LCI) / 3.92`.
- **`as_loftool`**: Placeholder column, always 0.0 (default). Not used in either Stan model.
- **ClinVar deduplication**: ClinVar evidence may have multiple rows per variant-disease (different submissions). Deduplication happens in Step 3 `process_gene()`, not Step 1.
- **as_clinicalSignificance**: Sourced from OT ClinVar evidence in Step 1 (preferred) or from VEP `--check_existing` in Step 2. Step 3 merges both, preferring the OT value when > 0.
- **FoldX transform**: Per-gene (not global) — computed in `process_gene()` using `1 - (inverted / max_inverted)` where inverted = `max_score - raw_score`.
- **Spark scratch dir**: Stan writes to `SPARK_LOCAL_DIRS` (SSD-backed), not `/tmp` (RAM-backed tmpfs), to avoid EPERM errors for large gene-disease pairs.
