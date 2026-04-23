"""VIDRA Data Preparation Pipeline (PySpark on Dataproc Serverless).

Replaces the original Nextflow + BigQuery pipeline with a single PySpark job
that reads all data from GCS and writes analysis-ready parquet partitioned by gene.

Pipeline steps:
  0.  Load gene mapping (Ensembl ID <-> Gene Symbol) from OT Platform 24.03 targets parquet,
      build obsolete-to-current EFO remapping from OT Platform 24.03 diseases parquet,
      and load study->EFO phenotype map
  1.  Load AZ EFO phenotype map from curated TSV (az_phewas_to_efo.tsv)
  2.  AZ PheWAS rare variants: allelic model, p<5e-8, gene symbol->ensembl join,
      EFO mapping, SE from log(OR CI). GsourceLab=1, GqtlLab=2.
  3.  ClinVar rare variants: datasourceId=eva, confidence filter, disease mapping.
      GsourceLab=2, GqtlLab=2.
  4.  Common variants: coloc H4>0.7, blood tissue filter, GWAS+QTL join,
      study->EFO mapping, dedup by best GWAS then QTL pval. GsourceLab=0.
  5.  Coding GWAS: OT variants (coding consequences), GWAS pval<=5e-8,
      anti-join to exclude coloc variants, dedup. GsourceLab=3, GqtlLab=2.
  6.  Burden tests: raw AZ PheWAS collapsing model data (binary + quantitative),
      gene symbol->ensembl join, EFO mapping, dedup by lowest pValue.
  7.  Union all 4 variant sources, left-join burden on gene+disease.
  8.  Fill missing values with defaults, write partitioned parquet.

Data sources (all read from gs://<bucket>/raw_data/):
  - AZ PheWAS CSV:    raw_data/2022_03_07_azphewas-com-450k-exwas-binary.csv.bz2
  - AZ EFO curation:  raw_data/az_phewas_to_efo.tsv (curated phenotype-to-EFO mapping)
  - Gene mapping:     raw_data/open_targets/targets/ (OT Platform 24.03 target index)
  - Disease mapping:  raw_data/open_targets/diseases/ (OT Platform 24.03 disease index, obsoleteTerms)
  - OT coloc:         raw_data/open_targets/coloc/       (v2d_coloc parquet)
  - OT GWAS:          raw_data/open_targets/gwas/        (sa_gwas parquet)
  - OT ClinVar:       raw_data/open_targets/clinvar/     (evidence, sourceId=eva)
  - AZ PheWAS burden:  raw_data/az_phewas_burden/{binary,quantitative}/ (collapsing model parquet)
  - OT variants:      raw_data/open_targets/variants/    (variant-index parquet)
  - OT studies:       raw_data/open_targets/studies/     (study-index parquet)

  See scripts/copy_opentargets_data.sh for how to copy OT data into the bucket.

Output:
  gs://<bucket>/vidra_analysis_ready/   (parquet, partitioned by as_gene)
  Consumed by: scripts/pyspark_scripts/run_bayesian_analysis.py

Output columns:
  variant, as_gene, as_disease, GsourceLab, GqtlLab,
  yc, ycse, xc, xcse, bO, bOse, as_clinicalSignificance

  Note: VEP-derived annotations (blosum62, sift, polyphen, cadd, revel, etc.)
  are NOT included here — they come from the Step 2 annotation parquet and are
  joined in run_bayesian_analysis.py. Only as_clinicalSignificance is carried
  from Step 1 because ClinVar's OT evidence provides a higher-quality value
  than VEP's --check_existing for that column.

Intentional improvements over the original Nextflow pipeline:
  - AZ SE: uses (log(UCI) - log(LCI)) / 3.92 instead of sqrt(50) * (UCI-LCI) / 3.92
  - Burden: uses raw AZ PheWAS collapsing model data (p<0.1) instead of OT Platform
    gene_burden evidence (pre-filtered to ~p<5e-8). Binary traits use log(OR) scale,
    quantitative traits use beta directly. SE = (UCI - LCI) / 3.92 on the appropriate
    scale. This rescues ~35k additional burden tests over the OT pre-filtered data.
  - Coding GWAS: includes 'inframe_deletion' (missing from original's gene-filtered query)


Known limitations:
  - VEP annotations: using GERP as the conservation score instead of conservScore from ProtVar.
  - QTL effects for common variants come from the coloc table
    (left_var_right_study_beta/se) rather than a separate sa_molecular_trait join.
    These should be near-identical since coloc implies shared causal variants.

Submit to Dataproc Serverless:
  # Full production run:
  gcloud dataproc batches submit pyspark \\
    scripts/pyspark_scripts/prepare_analysis_input.py \\
    --project=open-targets-genetics-dev \\
    --region=europe-west1 \\
    --deps-bucket=gs://vidra-2-0 \\
    --container-image=europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \\
    --properties=spark.sql.execution.arrow.pyspark.enabled=true \\
    -- --bucket_name=vidra-2-0

  # Test mode (200 random genes):
  gcloud dataproc batches submit pyspark \\
    scripts/pyspark_scripts/prepare_analysis_input.py \\
    --project=open-targets-genetics-dev \\
    --region=europe-west1 \\
    --deps-bucket=gs://vidra-2-0 \\
    --container-image=europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \\
    --properties=spark.sql.execution.arrow.pyspark.enabled=true \\
    -- --bucket_name=vidra-2-0 --test_mode --test_genes 200

  # Check job logs:
  BATCH_ID=<id from submit output>
  OUTPUT_URI=$(gcloud dataproc batches describe $BATCH_ID \\
    --project=open-targets-genetics-dev --region=europe-west1 \\
    --format="value(runtimeInfo.outputUri)")
  gsutil cat "${OUTPUT_URI}.000000000"
"""

from pyspark.sql import SparkSession, Window
from pyspark.sql import functions as F
from pyspark.sql.types import StructType, StructField, StringType, DoubleType


# ============================================================================
# Constants
# ============================================================================
COLOC_H4_THRESHOLD = 0.7
AZ_PVAL_THRESHOLD = 5e-8

# ClinVar clinical significance ordinal encoding.
# Matches original pre_processing_VIDRA_per_gene_pheno.py (lines 458-468).
# Order: 0 = not provided / benign → 16 = pathogenic.  MinMaxScale to [0, 1].
CLINICAL_SIG_ORDER = [
    'not provided', 'association not found', 'other', 'benign', 'likely benign',
    'low penetrance', 'confers sensitivity', 'uncertain risk allele',
    'drug response',
    'uncertain significance', 'association', 'affects',
    'likely risk allele', 'risk factor',
    'established risk allele', 'likely pathogenic', 'pathogenic',
]
_N_CLIN_CATS = len(CLINICAL_SIG_ORDER) - 1   # 16, denominator for MinMaxScale
_CLIN_SIG_LOOKUP = {cat: i / _N_CLIN_CATS for i, cat in enumerate(CLINICAL_SIG_ORDER)}
# OT ClinVar evidence emits "conflicting interpretations of pathogenicity" as
# a composite clinicalSignificances value (~38k rows). It isn't on the 17-term
# ordinal above, so without this alias it would fall through to 0.0 (benign),
# which is the opposite of what conflicting evidence implies. Score it as
# equivalent to "uncertain significance" — genuinely undetermined.
_CLIN_SIG_LOOKUP['conflicting interpretations of pathogenicity'] = \
    _CLIN_SIG_LOOKUP['uncertain significance']

# Blood-related tissues to keep for common QTL variants
# From original: pre_processing_VIDRA_per_gene_pheno.py
BLOOD_TISSUES = [
    'UBERON_0000178', 'BLOOD', 'UBERON_0001969',
    'TREG_NAIVE', 'TREG_MEMORY', 'TH2_MEMORY', 'TH1_MEMORY', 'TH17_MEMORY',
    'TH1-17_MEMORY', 'TFH_MEMORY', 'NK-CELL_NAIVE',
    'MONOCYTE_NAIVE', 'MONOCYTE_CD16_NAIVE',
    'CD8_T-CELL_NAIVE', 'CD8_T-CELL_ANTI-CD3-CD28',
    'CD4_T-CELL_NAIVE', 'CD4_T-CELL_ANTI-CD3-CD28',
    'B-CELL_NAIVE',
    'MONOCYTE_R848', 'MONOCYTE_PAM3CSK4', 'MONOCYTE_LPS',
    'MONOCYTE_IAV', 'MACROPHAGE_SALMONELLA', 'MACROPHAGE_NAIVE',
    'MACROPHAGE_LISTERIA',
    'NEUTROPHIL_CD16', 'T-CELL_CD8', 'T-CELL_CD4', 'NEUTROPHIL',
    'MONOCYTE', 'MACROPHAGE_IFNG', 'MACROPHAGE_IFNG+SALMONELLA',
    'T-CELL', 'MONOCYTE_LPS24', 'MONOCYTE_LPS2', 'MONOCYTE_IFN24',
    'B-CELL_CD19', 'PLATELET', 'NEUTROPHIL_CD15', 'MONOCYTE_CD14',
]

# Coding consequences for GWAS coding variants
CODING_CONSEQUENCES = [
    'inframe_insertion', 'frameshift_variant', 'stop_gained',
    'splice_donor_variant', 'coding_sequence_variant', 'stop_lost',
    'stop_retained_variant', 'missense_variant',
    'incomplete_terminal_codon_variant', 'protein_altering_variant',
    'start_lost', 'synonymous_variant', 'splice_acceptor_variant',
    'inframe_deletion',
]

# Burden test collapsing models (raw AZ PheWAS data, pre-filtered to p<0.1)
BURDEN_COLLAPSING_MODELS = ['ptv', 'ptv5pcnt', 'ptvraredmg']

# AZ CSV Schema (bzip2 compressed)
AZ_SCHEMA = StructType([
    StructField("Variant", StringType(), True),
    StructField("Variant type", StringType(), True),
    StructField("Phenotype", StringType(), True),
    StructField("Category", StringType(), True),
    StructField("Model", StringType(), True),
    StructField("Consequence type", StringType(), True),
    StructField("Gene", StringType(), True),
    StructField("Transcript", StringType(), True),
    StructField("cDNA change", StringType(), True),
    StructField("Amino acid change", StringType(), True),
    StructField("Exon rank", StringType(), True),
    StructField("No. cases", DoubleType(), True),
    StructField("No. AA cases", DoubleType(), True),
    StructField("No. AB cases", DoubleType(), True),
    StructField("No. BB cases", DoubleType(), True),
    StructField("Case MAF", DoubleType(), True),
    StructField("% AB or BB cases", DoubleType(), True),
    StructField("% BB cases", DoubleType(), True),
    StructField("No. controls", DoubleType(), True),
    StructField("No. AA controls", DoubleType(), True),
    StructField("No. AB controls", DoubleType(), True),
    StructField("No. BB controls", DoubleType(), True),
    StructField("Control MAF", DoubleType(), True),
    StructField("% AB or BB controls", DoubleType(), True),
    StructField("% BB controls", DoubleType(), True),
    StructField("p-value", DoubleType(), True),
    StructField("Odds ratio", DoubleType(), True),
    StructField("Odds ratio LCI", DoubleType(), True),
    StructField("Odds ratio UCI", DoubleType(), True),
])

# Default fill values for missing effect-size data.
# Real VEP annotations are joined in run_bayesian_analysis.py from the
# Step 2 annotation parquet — only as_clinicalSignificance is carried
# from Step 1 (ClinVar OT evidence value, merged with VEP in Step 3).
DEFAULT_FILL = {
    'as_clinicalSignificance': 0.0,
    'bO': 0.0, 'bOse': 2.0,
    'ycse': 0.14, 'yc': 0.0,
    'xc': 0.0, 'xcse': 0.1,
}

# Final output columns
OUTPUT_COLS = [
    'variant', 'as_gene', 'as_disease',
    'GsourceLab', 'GqtlLab',
    'yc', 'ycse', 'xc', 'xcse', 'bO', 'bOse',
    'as_clinicalSignificance',
]


def main(args):
    spark = SparkSession.builder \
        .appName("VIDRA_Prepare_Analysis_Input") \
        .getOrCreate()

    BUCKET = args.bucket_name
    RAW_ROOT = f"gs://{BUCKET}/raw_data"
    OT_ROOT = f"{RAW_ROOT}/open_targets"
    OUTPUT_ROOT = f"gs://{BUCKET}/vidra_analysis_ready{args.output_suffix}"
    TEST_MODE = args.test_mode
    TEST_GENES = args.test_genes

    print(f"--- Starting Analysis Prep (test_mode={TEST_MODE}) ---")

    # =========================================================================
    # 0. REFERENCE DATA
    # =========================================================================
    print("Loading Gene Mapping...")
    gene_map_df = spark.read.parquet(f"{OT_ROOT}/targets/") \
        .select(
            F.col("id").alias("ensembl_id"),
            F.col("approvedSymbol").alias("gene_symbol")
        ) \
        .filter(F.col("ensembl_id").isNotNull())

    if TEST_MODE:
        # Random sample of N genes (avoid first-N bias toward mitochondrial genes)
        gene_map_df = gene_map_df.orderBy(F.rand(seed=42)).limit(TEST_GENES)
        gene_map_df.cache()
        print(f"TEST MODE: Limited to {gene_map_df.count()} genes")

    gene_ensembl_list = gene_map_df.select("ensembl_id")

    # =========================================================================
    # 0b. OBSOLETE EFO ID REMAPPING (from OT Platform diseases index)
    # =========================================================================
    # The OT diseases parquet has an `obsoleteTerms` array listing old EFO IDs
    # that have been replaced by each current disease ID. Exploding this gives a
    # lookup table (old_efo -> current_efo) used to normalise disease IDs across
    # all sources before any cross-source joins.
    # This replaces the legacy convert_old_EFOids_toCurrentOnes.csv from the
    # original Nextflow pipeline.
    print("Loading Disease EFO Obsolete-to-Current Mapping...")
    diseases_df = spark.read.parquet(f"{OT_ROOT}/diseases/")
    efo_remap = diseases_df \
        .select(
            F.explode("obsoleteTerms").alias("old_efo"),
            F.col("id").alias("current_efo"),
        ) \
        .dropDuplicates(["old_efo"])
    print(f"EFO Remap entries: {efo_remap.count()}")

    def apply_efo_remap(df, disease_col):
        """Replace obsolete EFO IDs in disease_col with their current equivalents."""
        return df \
            .join(efo_remap, df[disease_col] == efo_remap["old_efo"], "left") \
            .withColumn(disease_col, F.coalesce(F.col("current_efo"), F.col(disease_col))) \
            .drop("old_efo", "current_efo")

    # =========================================================================
    # 0c. STUDY -> PHENOTYPE MAPPING (from OT Genetics study index)
    # =========================================================================
    # Maps GWAS study_id to EFO disease IDs (trait_efos is an array)
    print("Loading Study -> Phenotype Mapping...")
    studies_raw = spark.read.parquet(f"{OT_ROOT}/studies/")
    # Explode trait_efos array to get one row per study_id + EFO
    phenoTab = studies_raw \
        .select("study_id", F.explode("trait_efos").alias("diseaseId")) \
        .dropDuplicates()
    # Remap obsolete EFOs from the 22.09 Genetics study index to current IDs
    phenoTab = apply_efo_remap(phenoTab, "diseaseId")
    print(f"Study-Phenotype Map entries: {phenoTab.count()}")

    # =========================================================================
    # 1. EFO PHENOTYPE MAPPING (from curated TSV)
    # =========================================================================
    # Curated mapping from phenotype strings to EFO disease IDs. Used for both
    # AZ PheWAS rare variants (Section 2) and AZ PheWAS burden (Section 6).
    #
    # We deliberately do NOT filter to STUDY contains "AstraZeneca" here: AZ
    # phenotype names are sometimes curated under other study labels (e.g.
    # "LDL direct" is curated under "Genebass"), so filtering would drop valid
    # mappings for both rare variants and burden. The previous code filtered
    # here for rare variants only, causing an asymmetry where burden had fuller
    # coverage than rare variants for the same phenotype strings.
    print("Loading EFO Phenotype Map...")
    efo_curation = spark.read \
        .option("header", "true") \
        .option("sep", "\t") \
        .csv(f"{RAW_ROOT}/az_phewas_to_efo.tsv")

    efo_map = efo_curation \
        .withColumn("diseaseId",
                    F.element_at(F.split(F.col("SEMANTIC_TAG"), "/"), -1)) \
        .select(
            F.col("PROPERTY_VALUE").alias("diseaseFromSource"),
            F.col("diseaseId"),
        ) \
        .dropDuplicates() \
        .filter(F.col("diseaseId").isNotNull())

    print(f"EFO Map entries: {efo_map.count()}")

    # =========================================================================
    # 2. AZ PHEWAS RARE VARIANTS
    # =========================================================================
    print("Loading AZ PheWAS Data...")
    az_raw_path = f"{RAW_ROOT}/2022_03_07_azphewas-com-450k-exwas-binary.csv.bz2"

    az_raw = spark.read \
        .option("header", "true") \
        .schema(AZ_SCHEMA) \
        .csv(az_raw_path)

    # Filter 1: Allelic model only (from original pipeline)
    az_filtered = az_raw.filter(F.col("Model") == "allelic")

    # Filter 2: Significant p-values
    az_filtered = az_filtered.filter(F.col("p-value") < AZ_PVAL_THRESHOLD)

    # Clean gene names: remove quotes (AZ artefact)
    az_filtered = az_filtered.withColumn(
        "Gene", F.regexp_replace(F.col("Gene"), "'", "")
    )

    # Filter 3: Keep only genes in our gene list
    # Join directly with gene_map_df to get both filter AND ensembl_id in one step
    az_mapped = az_filtered.join(
        gene_map_df,
        az_filtered.Gene == gene_map_df.gene_symbol,
        "inner"
    ).drop(gene_map_df.gene_symbol)  # Drop duplicate column from right side

    if TEST_MODE:
        print(f"AZ Filtered Count (test): {az_mapped.count()}")
    else:
        print(f"AZ Filtered Count: {az_mapped.count()}")

    # Map phenotypes to EFO IDs
    az_with_efo = az_mapped.join(
        efo_map,
        az_mapped.Phenotype == efo_map.diseaseFromSource,
        "left"
    ).drop(efo_map.diseaseFromSource)  # Drop duplicate column

    # Use EFO ID if available; drop rows without EFO mapping.
    # Raw phenotype names would never match ClinVar/burden (which use EFO IDs).
    az_with_efo = az_with_efo.filter(
        F.col("diseaseId").isNotNull()
    ).withColumn("as_disease", F.col("diseaseId"))
    # Remap any obsolete EFO IDs from the curated TSV to current ones
    az_with_efo = apply_efo_remap(az_with_efo, "as_disease")

    # Harmonise variant IDs: replace hyphens/spaces with underscores
    az_with_efo = az_with_efo.withColumn(
        "variant_clean",
        F.regexp_replace(F.regexp_replace(F.col("Variant"), "-", "_"), " ", "_")
    )

    # Calculate Beta and SE (corrected formula)
    # beta = log(OR)
    # SE = (log(UCI) - log(LCI)) / 3.92
    az_with_efo = az_with_efo \
        .withColumn("yc", F.log(F.col("Odds ratio"))) \
        .withColumn("ycse",
                    (F.log(F.col("Odds ratio UCI")) - F.log(F.col("Odds ratio LCI"))) / 3.92)

    # Select standardised output columns
    az_ready = az_with_efo.select(
        F.col("variant_clean").alias("variant"),
        F.col("ensembl_id").alias("as_gene"),
        F.col("as_disease"),
        F.lit(1).cast("int").alias("GsourceLab"),   # 1 = AZ_PheWAS
        F.lit(2).cast("int").alias("GqtlLab"),      # 2 = no_qtl
        F.col("yc"),
        F.col("ycse"),
        F.lit(0.0).alias("xc"),
        F.lit(0.1).alias("xcse"),
        F.lit(0.0).alias("as_clinicalSignificance"),
    ).dropDuplicates()  # matches original drop_duplicates()

    print(f"AZ Ready Count: {az_ready.count()}")

    # =========================================================================
    # 3. CLINVAR RARE VARIANTS
    # =========================================================================
    print("Loading ClinVar Data...")
    clinvar_raw = spark.read.parquet(f"{OT_ROOT}/clinvar/")

    # Filter: datasourceId=eva, target in gene list, adequate confidence
    clinvar_filtered = clinvar_raw \
        .filter(F.col("datasourceId") == "eva") \
        .filter(F.col("variantId").isNotNull()) \
        .filter(~F.col("confidence").isin(
            # "criteria provided, single submitter",
            "criteria provided, conflicting interpretations",
            "no assertion criteria provided",
            "no assertion provided"
        )) \
        .join(gene_ensembl_list,
              clinvar_raw.targetId == gene_ensembl_list.ensembl_id,
              "inner") \
        .drop(gene_ensembl_list.ensembl_id)  # Drop duplicate column

    # Extract clinical significance from OT evidence (same source as original).
    # clinicalSignificances is an array column; take the first element.
    clinvar_with_sig = clinvar_filtered.withColumn(
        "_clin_sig_raw",
        F.lower(F.trim(F.col("clinicalSignificances").getItem(0)))
    )

    # Build a Spark mapping expression for the 17-category ordinal → [0, 1]
    _clin_sig_mapping = F.create_map(
        *[item for pair in [
            (F.lit(cat), F.lit(score))
            for cat, score in _CLIN_SIG_LOOKUP.items()
        ] for item in pair]
    )

    clinvar_with_sig = clinvar_with_sig.withColumn(
        "_as_clinicalSignificance",
        F.coalesce(_clin_sig_mapping[F.col("_clin_sig_raw")], F.lit(0.0))
    )

    clinvar_df = clinvar_with_sig.select(
        F.col("variantId").alias("variant"),
        F.col("targetId").alias("as_gene"),
        F.col("diseaseId").alias("as_disease"),
        F.lit(2).cast("int").alias("GsourceLab"),   # 2 = ClinVar
        F.lit(2).cast("int").alias("GqtlLab"),      # 2 = no_qtl
        # ClinVar has no GWAS effect and no QTL effect — use same defaults as
        # original pre_processing_VIDRA_per_gene_pheno.py (beta_gwas=0, se_gwas=0.14,
        # beta_qtl=0, se_qtl=0.1). These values are not used by the Stan ClinVar
        # branch (which only uses disease_prior / protein_prior), but aligning with
        # the original avoids any difference via the vectorised xc~normal(xcest, xcse)
        # measurement model that runs over all N rows when source-0 variants are present.
        F.lit(0.0).alias("yc"),
        F.lit(DEFAULT_FILL['ycse']).alias("ycse"),  # 0.14
        F.lit(0.0).alias("xc"),
        F.lit(0.1).alias("xcse"),
        # ClinVar clinical significance from OT evidence (NOT VEP).
        # This is the only annotation carried from Step 1 — it gets merged
        # with VEP's as_clinicalSignificance in Step 3, preferring this value.
        F.col("_as_clinicalSignificance").alias("as_clinicalSignificance"),
    )

    # Drop rows without a mapped disease
    clinvar_df = clinvar_df.filter(F.col("as_disease").isNotNull())
    # Remap any obsolete EFO IDs from the 24.03 ClinVar evidence to current ones
    clinvar_df = apply_efo_remap(clinvar_df, "as_disease")

    print(f"ClinVar Ready Count: {clinvar_df.count()}")

    # =========================================================================
    # 4. COMMON VARIANTS (Coloc + GWAS + QTL)
    # =========================================================================
    print("Loading Colocalisation Data...")
    coloc_raw = spark.read.parquet(f"{OT_ROOT}/coloc/")

    # Filter coloc: H4 > threshold, exclude gwas-gwas and sqtl
    coloc_filtered = coloc_raw \
        .filter(F.col("coloc_h4") > COLOC_H4_THRESHOLD) \
        .filter(F.col("right_type") != "gwas") \
        .filter(F.col("right_type") != "sqtl") \
        .filter(F.col("left_type") == "gwas") \
        .join(gene_ensembl_list,
              coloc_raw.right_gene_id == gene_ensembl_list.ensembl_id,
              "inner") \
        .drop(gene_ensembl_list.ensembl_id)  # Drop duplicate column

    # Blood tissue filter (from original pipeline)
    coloc_filtered = coloc_filtered.filter(
        F.col("right_bio_feature").isin(BLOOD_TISSUES)
    )

    # Drop coloc rows where the GWAS lead variant's effect in the QTL study
    # wasn't propagated by OT (NULL left_var_right_study_beta/se). Without this,
    # the dedup window's asc() ordering is NULLS FIRST, so a NULL-beta row
    # silently beats a valid sibling for the same (variant, gene, disease) and
    # xc ends up as a no-effect sentinel — fabricating a no-QTL data point for
    # a coloc hit. ~6% of blood-filtered coloc rows; see null_qtl_beta_diagnostic
    # notebook for impact analysis.
    coloc_filtered = coloc_filtered.filter(
        F.col("left_var_right_study_beta").isNotNull()
        & F.col("left_var_right_study_se").isNotNull()
    )

    print(f"Coloc Filtered Count: {coloc_filtered.count()}")

    # Load GWAS summary stats
    print("Loading GWAS Summary Stats...")
    gwas_raw = spark.read.parquet(f"{OT_ROOT}/gwas/")

    # Join GWAS onto coloc (left_study+left_chrom+left_pos+left_ref+left_alt)
    common_with_gwas = coloc_filtered.join(
        gwas_raw,
        (coloc_filtered.left_study == gwas_raw.study_id) &
        (coloc_filtered.left_chrom == gwas_raw.chrom) &
        (coloc_filtered.left_pos == gwas_raw.pos) &
        (coloc_filtered.left_ref == gwas_raw.ref) &
        (coloc_filtered.left_alt == gwas_raw.alt),
        "inner"
    ).drop(gwas_raw.study_id).drop(gwas_raw.chrom).drop(gwas_raw.pos) \
     .drop(gwas_raw.ref).drop(gwas_raw.alt)

    # Map GWAS study_id to EFO disease using the study index
    common_with_gwas = common_with_gwas.join(
        phenoTab,
        common_with_gwas.left_study == phenoTab.study_id,
        "inner"
    ).drop(phenoTab.study_id)

    # Construct variant ID
    common_with_gwas = common_with_gwas.withColumn(
        "variantId",
        F.concat_ws("_",
                     F.col("left_chrom"),
                     F.col("left_pos").cast("string"),
                     F.col("left_ref"),
                     F.col("left_alt"))
    )

    # QTL effects come from the coloc table itself:
    # left_var_right_study_beta, left_var_right_study_se
    # These are the QTL effects of the GWAS lead variant in the QTL study

    # Deduplicate: per (variantId, right_gene_id, disease),
    # keep best GWAS pval, then best QTL pval.
    #
    # The original pipeline (pre_processing_VIDRA_per_gene_pheno.py) applies
    # two sequential filters: best_gwas (keeps ALL ties on min GWAS pval),
    # then best_qtl (keeps ALL ties on min QTL pval), then drop_duplicates.
    # This is equivalent to a single composite sort: primary key = GWAS pval,
    # secondary key = QTL pval, picking the first row per group.
    gwas_pval_col = F.col("pval")  # from GWAS table
    qtl_pval_col = F.col("left_var_right_study_pval")

    w_dedup = Window.partitionBy("variantId", "right_gene_id", "diseaseId") \
        .orderBy(gwas_pval_col.asc_nulls_last(), qtl_pval_col.asc_nulls_last())
    common_deduped = common_with_gwas \
        .withColumn("rn", F.row_number().over(w_dedup)) \
        .filter(F.col("rn") == 1) \
        .drop("rn")

    # Select standardised output
    common_ready = common_deduped.select(
        F.col("variantId").alias("variant"),
        F.col("right_gene_id").alias("as_gene"),
        F.col("diseaseId").alias("as_disease"),
        F.lit(0).cast("int").alias("GsourceLab"),  # 0 = common_QTL
        F.when(F.col("right_type") == "eqtl", 0)
         .otherwise(1).cast("int").alias("GqtlLab"),
        # GWAS effect
        F.col("beta").alias("yc"),
        F.col("se").alias("ycse"),
        # QTL effect (from coloc table; guaranteed non-null by filter above)
        F.col("left_var_right_study_beta").alias("xc"),
        F.col("left_var_right_study_se").alias("xcse"),
        F.lit(0.0).alias("as_clinicalSignificance"),
    )

    print(f"Common Variants Ready Count: {common_ready.count()}")

    # =========================================================================
    # 5. CODING GWAS VARIANTS (non-colocalising missense)
    # =========================================================================
    # From original: get_coding_GWAS_nonColoc.py
    # 1. Load OT variants table, filter to coding consequences with gene_id_prot_coding
    # 2. Join with GWAS stats on chr/pos/ref/alt (pval < 5e-8)
    # 3. Attach study -> EFO phenotype mapping
    # 4. Exclude variants already in coloc set
    # 5. GsourceLab=3, GqtlLab=2
    print("Loading Coding GWAS Variants...")
    variants_raw = spark.read.parquet(f"{OT_ROOT}/variants/")

    # Filter to coding variants with a protein-coding gene assignment
    coding_variants = variants_raw \
        .filter(F.col("gene_id_prot_coding").isNotNull()) \
        .filter(F.col("most_severe_consequence").isin(CODING_CONSEQUENCES))

    # Filter to our gene list
    coding_variants = coding_variants.join(
        gene_ensembl_list,
        coding_variants.gene_id_prot_coding == gene_ensembl_list.ensembl_id,
        "inner"
    ).drop(gene_ensembl_list.ensembl_id)

    print(f"Coding variants (in gene list): {coding_variants.count()}")

    # Join with GWAS summary stats on position (pval < 5e-8)
    gwas_sig = gwas_raw.filter(F.col("pval") <= AZ_PVAL_THRESHOLD)

    coding_gwas = coding_variants.join(
        gwas_sig,
        (coding_variants.chr_id == gwas_sig.chrom) &
        (coding_variants.position == gwas_sig.pos) &
        (coding_variants.ref_allele == gwas_sig.ref) &
        (coding_variants.alt_allele == gwas_sig.alt),
        "inner"
    )

    print(f"Coding GWAS matched: {coding_gwas.count()}")

    # Attach study -> EFO phenotype mapping
    coding_gwas = coding_gwas.join(
        phenoTab,
        coding_gwas.study_id == phenoTab.study_id,
        "inner"
    ).drop(phenoTab.study_id)

    # Construct variant ID
    coding_gwas = coding_gwas.withColumn(
        "variantId",
        F.concat_ws("_",
                     F.col("chr_id"),
                     F.col("position").cast("string"),
                     F.col("ref_allele"),
                     F.col("alt_allele"))
    )

    # Exclude variants already in the coloc set (anti-join on variantId)
    # Build a set of coloc variant IDs for exclusion
    coloc_variant_ids = common_ready.select(
        F.col("variant").alias("coloc_variant")
    ).distinct()

    coding_gwas = coding_gwas.join(
        coloc_variant_ids,
        coding_gwas.variantId == coloc_variant_ids.coloc_variant,
        "left_anti"
    )

    print(f"Coding GWAS after coloc exclusion: {coding_gwas.count()}")

    # Deduplicate: per (variantId, gene, disease) keep best GWAS pval
    w_cd_gwas = Window.partitionBy(
        "variantId", "gene_id_prot_coding", "diseaseId"
    ).orderBy(F.col("pval").asc())
    coding_gwas_deduped = coding_gwas \
        .withColumn("rn", F.row_number().over(w_cd_gwas)) \
        .filter(F.col("rn") == 1) \
        .drop("rn")

    # Select standardised output
    coding_gwas_ready = coding_gwas_deduped.select(
        F.col("variantId").alias("variant"),
        F.col("gene_id_prot_coding").alias("as_gene"),
        F.col("diseaseId").alias("as_disease"),
        F.lit(3).cast("int").alias("GsourceLab"),  # 3 = coding_GWAS
        F.lit(2).cast("int").alias("GqtlLab"),      # 2 = no_qtl
        # GWAS effect
        F.col("beta").alias("yc"),
        F.col("se").alias("ycse"),
        # No QTL for coding GWAS variants
        F.lit(0.0).alias("xc"),
        F.lit(0.1).alias("xcse"),
        F.lit(0.0).alias("as_clinicalSignificance"),
    )

    print(f"Coding GWAS Ready Count: {coding_gwas_ready.count()}")

    # =========================================================================
    # 6. BURDEN TESTS (raw AZ PheWAS collapsing model data)
    # =========================================================================
    # Uses the full raw AZ PheWAS collapsing model results (binary + quantitative)
    # instead of the OT Platform gene_burden evidence, which was pre-filtered to
    # ~p<5e-8 and yielded very few burden pairs. The raw data is pre-filtered by
    # the provider to p<0.1 and includes ~6.8M gene-disease pairs with EFO mapping.
    # Even weak burden tests (large SE) are more informative than the default
    # (bO=0.0, bOse=2.0) since they are at least data-informed.
    print("Loading Raw AZ PheWAS Burden Tests...")
    burden_raw = spark.read.option("mergeSchema", "true").parquet(
        f"{RAW_ROOT}/az_phewas_burden/azphewas-com-450k-phewas-binary/",
        f"{RAW_ROOT}/az_phewas_burden/azphewas-com-450k-phewas-quantitative/",
    )

    # Filter to relevant collapsing models
    burden_filtered = burden_raw.filter(
        F.col("CollapsingModel").isin(BURDEN_COLLAPSING_MODELS)
    )

    # Map gene symbols to Ensembl IDs and filter to our gene list
    burden_mapped = burden_filtered.join(
        gene_map_df,
        burden_filtered.Gene == gene_map_df.gene_symbol,
        "inner"
    ).drop(gene_map_df.gene_symbol)

    # Map phenotypes to EFO IDs using the shared efo_map from Section 1.
    # Phenotypes mapping to multiple EFOs (~4% of the curation file) produce
    # one burden row per EFO synonym — matching the rare-variant behaviour.
    burden_mapped = burden_mapped.join(
        efo_map,
        burden_mapped.Phenotype == efo_map.diseaseFromSource,
        "inner"
    ).drop(efo_map.diseaseFromSource)

    burden_mapped = burden_mapped.withColumnRenamed("diseaseId", "disease_key")

    # Compute burden effect size (bO) and SE (bOse).
    # Binary and quantitative parquets have different schemas:
    #   Binary:       BinOddsRatio, BinOddsRatioLCI, BinOddsRatioUCI (raw OR scale)
    #   Quantitative: beta, LCI, UCI (effect size scale)
    # Both need to produce log-scale estimates for the Stan model.
    burden_mapped = burden_mapped \
        .withColumn("burden_beta",
            F.when(F.col("Type") == "Quantitative", F.col("beta"))
             .otherwise(F.log(F.col("BinOddsRatio")))) \
        .withColumn("burden_se",
            F.when(F.col("Type") == "Quantitative",
                   (F.col("UCI") - F.col("LCI")) / 3.92)
             .otherwise(
                   (F.log(F.col("BinOddsRatioUCI")) - F.log(F.col("BinOddsRatioLCI"))) / 3.92))

    # Drop rows where effect size or SE could not be computed
    burden_mapped = burden_mapped.filter(
        F.col("burden_beta").isNotNull()
        & F.col("burden_se").isNotNull()
        & (F.col("burden_se") > 0)
    )

    # Keep best burden per gene-disease (lowest pValue)
    w_burden = Window.partitionBy("ensembl_id", "disease_key") \
        .orderBy(F.col("pValue").asc())
    burden_best = burden_mapped \
        .withColumn("rn", F.row_number().over(w_burden)) \
        .filter(F.col("rn") == 1) \
        .drop("rn")

    # Select burden columns for join
    burden_for_join = burden_best.select(
        F.col("ensembl_id").alias("b_gene"),
        F.col("disease_key").alias("b_disease"),
        F.col("burden_beta").alias("bO"),
        F.col("burden_se").alias("bOse"),
    )
    # Remap any obsolete EFO IDs from the curation file to current ones
    burden_for_join = apply_efo_remap(burden_for_join, "b_disease")

    # Dedup again after EFO remap: two different disease_keys may merge into the
    # same current EFO, creating duplicate (gene, disease) entries that would
    # multiply variant rows on the left join in Section 8.
    w_burden_post = Window.partitionBy("b_gene", "b_disease") \
        .orderBy(F.abs(F.col("bO")).desc())
    burden_for_join = burden_for_join \
        .withColumn("rn", F.row_number().over(w_burden_post)) \
        .filter(F.col("rn") == 1) \
        .drop("rn")

    print(f"Burden Tests Ready: {burden_for_join.count()}")

    # =========================================================================
    # 7. UNION ALL SOURCES
    # =========================================================================
    print("Union all sources...")

    # Union: AZ + ClinVar + Common + Coding GWAS
    all_variants = az_ready \
        .unionByName(clinvar_df, allowMissingColumns=True) \
        .unionByName(common_ready, allowMissingColumns=True) \
        .unionByName(coding_gwas_ready, allowMissingColumns=True)

    print(f"Total variants before burden join: {all_variants.count()}")

    # =========================================================================
    # 8. ATTACH BURDEN TESTS
    # =========================================================================
    print("Attaching burden tests...")
    # Left join burden on gene + disease
    all_with_burden = all_variants.join(
        burden_for_join,
        (all_variants.as_gene == burden_for_join.b_gene) &
        (all_variants.as_disease == burden_for_join.b_disease),
        "left"
    )

    # Use burden values if available, otherwise use defaults
    all_with_burden = all_with_burden \
        .withColumn("bO",
                    F.coalesce(F.col("bO"), F.lit(DEFAULT_FILL['bO']))) \
        .withColumn("bOse",
                    F.coalesce(F.col("bOse"), F.lit(DEFAULT_FILL['bOse']))) \
        .drop("b_gene", "b_disease")

    # =========================================================================
    # 9. FINAL CLEANING
    # =========================================================================
    print("Final cleaning...")

    # Fill NAs with defaults
    for col_name, default_val in DEFAULT_FILL.items():
        if col_name in [c.name for c in all_with_burden.schema]:
            all_with_burden = all_with_burden.fillna(
                {col_name: default_val}
            )

    # Remove rows missing gene or disease
    final_df = all_with_burden \
        .filter(F.col("as_gene").isNotNull()) \
        .filter(F.col("as_disease").isNotNull())

    # Select only the output columns
    final_df = final_df.select(*OUTPUT_COLS)

    # =========================================================================
    # 10. WRITE OUTPUT
    # =========================================================================
    # Repartition so each gene gets exactly ONE parquet file.
    # Without this, Spark's default 200 shuffle partitions scatter each gene's
    # data across ~10-30 small files, creating ~300K files total.  That makes
    # downstream reads (Step 3 Spark partition pruning) painfully
    # slow because every read must discover and open hundreds of thousands of
    # tiny files over GCS.
    n_genes = final_df.select("as_gene").distinct().count()
    print(f"Repartitioning to {n_genes} partitions (1 per gene)...")
    final_df = final_df.repartition(n_genes, "as_gene")

    print("Writing output...")
    final_df.write.mode("overwrite").partitionBy("as_gene").parquet(OUTPUT_ROOT)

    # Write a single-file manifest of unique variant IDs for annotate_variants.py
    manifest_path = f"gs://{args.bucket_name}/vidra_analysis_ready{args.output_suffix}_manifest"
    print("Writing unique variant manifest...")
    final_df.select("variant").distinct().coalesce(1).write.mode("overwrite").parquet(manifest_path)
    print(f"  Variant manifest: {manifest_path}")

    gene_count = n_genes
    total_count = final_df.count()
    print(f"Preparation Complete.")
    print(f"  Total variants: {total_count}")
    print(f"  Unique genes: {gene_count}")
    print(f"  Output: {OUTPUT_ROOT}")

    spark.stop()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="VIDRA Data Preparation Pipeline"
    )
    parser.add_argument("--bucket_name", required=True,
                        help="GCS bucket name (without gs:// prefix)")
    parser.add_argument("--test_mode", action="store_true",
                        help="Run in test mode with limited genes")
    parser.add_argument("--test_genes", type=int, default=50,
                        help="Number of genes to use in test mode (default: 50)")
    parser.add_argument("--output_suffix", type=str, default="",
                        help="Suffix for output directory name (e.g. '_dev' -> vidra_analysis_ready_dev)")
    args = parser.parse_args()
    main(args)
