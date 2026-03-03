"""VIDRA Bayesian Analysis — PySpark implementation.

Reads analysis-ready data (Hive-partitioned by as_gene, 1 parquet file per gene),
joins variant annotations (38 MB broadcast table), applies direction inversions
and normalisations, then runs per-(gene, disease) Bayesian hierarchical models
via CmdStanPy variational inference and writes posterior summaries to parquet.

Pipeline position:
  1. prepare_analysis_input.py → gs://<bucket>/vidra_analysis_ready/  (1 file per gene)
  2. annotate_variants_cli.py  → gs://<bucket>/variant_annotations/   (38 MB parquet)
  3. THIS SCRIPT               → gs://<bucket>/vidra_results/

When using --gene_list, Spark reads ONLY the requested gene partitions via
literal .isin() filter (guaranteed partition pruning — no 300-executor overhead).
The 38 MB annotation table is broadcast-joined automatically by Spark.

Stan models (compiled in Docker image):
  - VIDRA.stan: multi-variant hierarchical model with slope_random[1..5] per source
  - VIDRA_single_variant.stan: single-variant model with one slope parameter

Output columns:
  gene, as_disease, parameter, model, n_variants, source, qtl,
  mean, median, pct_1, pct_2_5, pct_5, pct_10, pct_25, pct_40, pct_50,
  pct_60, pct_75, pct_90, pct_95, pct_97_5, pct_99, pp_slope_pos, pp_slope_neg

Submit (production):
  gcloud dataproc batches submit pyspark scripts/pyspark_scripts/run_bayesian_analysis.py \\
    --project=open-targets-genetics-dev --region=europe-west1 \\
    --deps-bucket=gs://vidra-2-0 \\
    --container-image=europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra-spark:v1 \\
    --properties=spark.sql.execution.arrow.pyspark.enabled=true \\
    -- --bucket_name=vidra-2-0

Submit (test — specific genes only):
  Upload a text file of ENSG IDs (one per line) to GCS, then:
  ... same as above ... -- --bucket_name=vidra-2-0 --gene_list=gs://vidra-2-0/my_genes.txt

Submit (test — specific genes + random padding to N total):
  Combines --gene_list with --test_mode. Guarantees the listed genes are included,
  then fills up to --test_genes total with random genes from the full dataset.
  ... same as above ... -- --bucket_name=vidra-2-0 --gene_list=gs://vidra-2-0/my_genes.txt --test_mode --test_genes 50

Submit (test — random genes only):
  ... same as above ... -- --bucket_name=vidra-2-0 --test_mode --test_genes 50
"""
from pyspark.sql import SparkSession
from pyspark.sql.types import (
    StructType, StructField, StringType, DoubleType, IntegerType
)
import pandas as pd
import numpy as np
import shutil
import tempfile
import os
from cmdstanpy import CmdStanModel

# =============================================================================
# Constants
# =============================================================================

# Maps (GsourceLab, GqtlLab) pairs → slope_random index in VIDRA.stan
# slope_random[1]=eQTL(0,0), [2]=pQTL(0,1), [3]=codingGWAS(3,2),
# [4]=AZ(1,2), [5]=ClinVar(2,2)
COMBINATION_SLOPE = {1: [0, 0], 2: [0, 1], 3: [3, 2], 4: [1, 2], 5: [2, 2]}

# Default values for annotation columns (used when annotation lookup is incomplete).
# Convention: 0 = pathogenic/damaging, 1 = benign/tolerated for protein annotations.
# This matches the original pre_processing_VIDRA_per_gene_pheno.py after its
# 1-x inversions of revel, polyphen, cadd, conservation.  Stan model sigma
# values (.05-.3) were calibrated on this convention.
ANNOTATION_DEFAULTS = {
    'as_blosum62': 1.0, 'as_conservation': 1.0, 'as_sift': 1.0,
    'as_polyphen': 1.0, 'as_cadd': 1.0, 'as_alphamissense': 1.0,
    'as_primateai': 0.0, 'as_revel': 1.0, 'as_clinicalSignificance': 0.0,
    'as_loftool': 0.0, 'as_plddt': 1.0, 'as_consequence': 1.0,
    # Defaults are filled AFTER 1-x inversions, so they are in protein_prior
    # convention (0=damaging, 1=benign). Missing → benign assumption.
    'as_foldx': 1.0,  # default for missing FoldX (1.0 = neutral)
    'foldxDdq_raw': float('nan'),  # raw FoldX ΔΔG for per-gene transform
}

# Annotation columns from annotate_variants_cli.py output.
ANNOTATION_COLS = [
    'as_blosum62', 'as_conservation', 'as_sift', 'as_polyphen',
    'as_cadd', 'as_alphamissense', 'as_revel', 'as_primateai',
    'as_loftool', 'as_plddt', 'as_consequence',
    'as_clinicalSignificance', 'foldxDdq_raw',
    'most_severe_consequence',
]

# Consequence types treated as loss-of-function: hardcode protein
# annotations to 0 (maximally damaging) because missense-specific
# tools (REVEL, AlphaMissense, FoldX) produce default (benign) values
# for these variant types.
# Matches original pipeline: pre_processing_VIDRA_per_gene_pheno.py:449
LOF_CONSEQUENCES = {
    'stop_gained', 'start_lost',
}

# Annotation columns that receive per-gene mean imputation and LoF
# hardcoding to 0 (maximally damaging).  Matches the original pipeline's
# colz = ['revel','cadd_phred','conservation','sift_score','polyphen_score'].
# We additionally include AlphaMissense and FoldX (new annotations).
GENE_MEAN_IMPUTE_COLS = [
    'as_revel', 'as_alphamissense', 'as_sift', 'as_polyphen',
    'as_foldx', 'as_cadd', 'as_conservation',
]

# Columns that need 1-x inversion (raw high = damaging, but convention
# is 0 = damaging, 1 = benign for the protein_prior in Stan).
# as_plddt: high pLDDT = well-structured region = mutation more consequential → invert.
INVERSION_COLS = ['as_revel', 'as_polyphen', 'as_cadd',
                  'as_alphamissense', 'as_plddt']

# Posterior summary column names (Spark-safe — no special chars)
POSTERIOR_COLS = [
    'mean', 'median', 'pct_1', 'pct_2_5', 'pct_5', 'pct_10', 'pct_25',
    'pct_40', 'pct_50', 'pct_60', 'pct_75', 'pct_90', 'pct_95',
    'pct_97_5', 'pct_99', 'pp_slope_pos', 'pp_slope_neg'
]

# Metadata columns in output
META_COLS = ['gene', 'as_disease', 'parameter', 'model',
             'n_variants', 'source', 'qtl']

ALL_OUTPUT_COLS = META_COLS + POSTERIOR_COLS

# =============================================================================
# Global Model Cache (lazy init on worker)
# =============================================================================
MODELS = {}


def _get_scratch_dir():
    """Return the Spark scratch directory for Stan temp files.

    Dataproc Serverless provides fast SSD scratch space via SPARK_LOCAL_DIRS.
    The default /tmp is often a size-capped tmpfs (RAM-backed), which causes
    'Operation not permitted' (EPERM) when Stan writes large output files
    for gene-disease pairs with many variants.

    Also sets TMPDIR so that CmdStanPy's internal tempfile.mkstemp (used for
    the data JSON) goes to the same scratch disk, not /tmp.
    """
    scratch = os.environ.get('SPARK_LOCAL_DIRS', '/tmp').split(',')[0]
    os.environ['TMPDIR'] = scratch
    return scratch


def get_models():
    """Lazily load pre-compiled Stan models on the worker node.

    Note: cpp_options are intentionally omitted here so that CmdStanPy uses
    the models pre-compiled during the Docker image build (at /opt/vidra/stan_models/).
    Passing any cpp_options would trigger a recompilation attempt, which fails
    because the spark user has no write permission to that directory.
    Spark handles gene-level parallelism externally, so Stan threading is not needed.
    """
    if not MODELS:
        print("Initializing Stan Models on Worker...")
        path_multi = '/opt/vidra/stan_models/VIDRA.stan'
        path_single = '/opt/vidra/stan_models/VIDRA_single_variant.stan'
        MODELS['multi'] = CmdStanModel(stan_file=path_multi)
        MODELS['single'] = CmdStanModel(stan_file=path_single)
        print("Stan Models Initialized.")
    return MODELS


# =============================================================================
# Posterior helpers
# =============================================================================

def clean_posteriorForAs(serie):
    """Compute summary stats for a posterior distribution.

    Returns dict with Spark-safe column names matching POSTERIOR_COLS.
    """
    s = np.asarray(serie)
    return {
        'mean': float(np.mean(s)),
        'median': float(np.median(s)),
        'pct_1': float(np.percentile(s, 1)),
        'pct_2_5': float(np.percentile(s, 2.5)),
        'pct_5': float(np.percentile(s, 5)),
        'pct_10': float(np.percentile(s, 10)),
        'pct_25': float(np.percentile(s, 25)),
        'pct_40': float(np.percentile(s, 40)),
        'pct_50': float(np.percentile(s, 50)),
        'pct_60': float(np.percentile(s, 60)),
        'pct_75': float(np.percentile(s, 75)),
        'pct_90': float(np.percentile(s, 90)),
        'pct_95': float(np.percentile(s, 95)),
        'pct_97_5': float(np.percentile(s, 97.5)),
        'pct_99': float(np.percentile(s, 99)),
        'pp_slope_pos': float((s > 0).mean()),
        'pp_slope_neg': float((s < 0).mean()),
    }


def _get_slope_key(combination_slope, value_to_find):
    """Return the dict key whose value matches value_to_find, or None."""
    for key, value in combination_slope.items():
        if value == value_to_find:
            return key
    return None


def _safe_float(val, default=0.0):
    """Convert to float, replacing NaN/None/inf with default."""
    try:
        f = float(val)
        return default if not np.isfinite(f) else f
    except (TypeError, ValueError):
        return default


def _fill_annotation(df, col, default=0.0):
    """Ensure column exists in df and fill NaN with default."""
    if col not in df.columns:
        df[col] = default
    else:
        df[col] = df[col].fillna(default)
    return df


# =============================================================================
# Stan model runners
# =============================================================================

def AS_singleVars(df, gene, h1=0.1):
    """Run single-variant Stan model (VIDRA_single_variant.stan).

    Args:
        df: pandas DataFrame with one row per variant (typically 1 row)
        gene: gene identifier string
        h1: hyperparameter (default 0.1)

    Returns:
        pandas DataFrame with one row of posterior summaries, or None on failure.
    """
    model = get_models()['single']
    if len(df) == 0:
        return None

    row = df.iloc[0]
    df_dict = {
        'h1': float(h1),
        'N': len(df),
        'numG1': _safe_float(row['GsourceLab']),
        'numG2': _safe_float(row['GqtlLab']),
        'xc': _safe_float(row['xc'], 0.0),
        'xcse': max(_safe_float(row['xcse'], 0.1), 1e-6),
        'yOR': _safe_float(row['yc']),
        'yORse': max(_safe_float(row['ycse'], 0.14), 1e-6),
        'as_blosum62': _safe_float(row.get('as_blosum62', 1.0)),
        'as_conservation': _safe_float(row.get('as_conservation', 1.0)),
        'as_sift': _safe_float(row.get('as_sift', 1.0)),
        'as_polyphen': _safe_float(row.get('as_polyphen', 1.0)),
        'as_clinicalSignificance': _safe_float(
            row.get('as_clinicalSignificance', 0.0)),
        'as_cadd': _safe_float(row.get('as_cadd', 1.0)),
        'as_alphamissense': _safe_float(row.get('as_alphamissense', 1.0)),
        'as_primateai': _safe_float(row.get('as_primateai', 0.0)),
    }

    # Use Spark's scratch SSD for Stan I/O — /tmp is a size-capped tmpfs
    # that causes EPERM for large gene-disease pairs.
    scratch = _get_scratch_dir()
    with tempfile.TemporaryDirectory(prefix='stan_sv_', dir=scratch) as tmpdir:
        fit = model.variational(
            data=df_dict, seed=412, algorithm='fullrank',
            grad_samples=20, draws=1000,
            require_converged=False, show_console=True, refresh=1000,
            output_dir=tmpdir
        )
        slope_posteriors = fit.stan_variable('slope', mean=False)
    res = clean_posteriorForAs(slope_posteriors)
    res['gene'] = gene
    res['as_disease'] = str(df['as_disease'].iloc[0])
    res['n_variants'] = int(df['variant'].nunique())
    res['source'] = str(int(df_dict['numG1']))
    res['qtl'] = str(int(df_dict['numG2']))
    res['model'] = 'single_variant'
    res['parameter'] = 'slope'
    return pd.DataFrame([res])


def AS_multiVars(df, gene):
    """Run multi-variant Stan model (VIDRA.stan).

    Faithfully replicates the original meta-slope extraction logic:
    1. Identify which (GsourceLab, GqtlLab) combinations are present
    2. Map them to slope_random indices via COMBINATION_SLOPE
    3. Pool posterior samples from observed slope_random components → meta_slope
    4. Return per-parameter posterior summaries + meta_slope row

    Args:
        df: pandas DataFrame with multiple variant rows for one disease
        gene: gene identifier string

    Returns:
        pandas DataFrame with one row per parameter, or None on failure.
    """
    model = get_models()['multi']

    nu = max(int(df['variant'].nunique()) - 1, 1)

    # Build Stan data dict — all list columns for N observations
    df_dict = {
        'nu': nu,
        'N': len(df),
        'numG1': df['GsourceLab'].astype(float).tolist(),
        'numG2': df['GqtlLab'].astype(float).tolist(),
        'xc': df['xc'].astype(float).replace([np.inf, -np.inf], np.nan).fillna(1.0).tolist(),
        'xcse': df['xcse'].astype(float).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(lower=1e-6).tolist(),
        'yOR': df['yc'].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0).tolist(),
        'yORse': df['ycse'].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.14).clip(lower=1e-6).tolist(),
        'bO': df['bO'].astype(float).tolist(),
        'bOse': df['bOse'].astype(float).tolist(),
        # Annotations used by VIDRA.stan protein_prior:
        'as_blosum62': df['as_blosum62'].astype(float).tolist(),
        'as_conservation': df['as_conservation'].astype(float).tolist(),
        'as_sift': df['as_sift'].astype(float).tolist(),
        'as_polyphen': df['as_polyphen'].astype(float).tolist(),
        'as_cadd': df['as_cadd'].astype(float).tolist(),
        'as_alphamissense': df['as_alphamissense'].astype(float).tolist(),
        'as_revel': df['as_revel'].astype(float).tolist(),
        # Annotations from VEP + FoldX (via annotate_variants_cli.py):
        'as_foldx': df['as_foldx'].astype(float).tolist(),
        'as_consequence': df['as_consequence'].astype(float).tolist(),
        # Disease priors:
        'as_clinicalSignificance': df['as_clinicalSignificance'].astype(
            float).tolist(),
        'as_primateai': df['as_primateai'].astype(float).tolist(),
    }

    # Use Spark's scratch SSD for Stan I/O — /tmp is a size-capped tmpfs
    # that causes EPERM for large gene-disease pairs.
    scratch = _get_scratch_dir()
    stan_tmpdir = tempfile.mkdtemp(prefix='stan_mv_', dir=scratch)

    # Try fullrank first; fall back to meanfield if fullrank fails.
    # Fullrank estimates a full covariance matrix over 4N+8 parameters,
    # which becomes ill-conditioned for large N (many variants).
    # Meanfield uses a diagonal covariance and is much more stable.
    algorithm = 'fullrank'
    try:
        fit = model.variational(
            data=df_dict, seed=412, algorithm=algorithm,
            grad_samples=20, draws=1000,
            require_converged=False, show_console=True, refresh=1000,
            output_dir=stan_tmpdir
        )
    except Exception as e_fullrank:
        print(f"  fullrank failed for {gene} (N={len(df)}): {e_fullrank}; "
              f"retrying with meanfield...")
        # Clean up failed output and retry
        shutil.rmtree(stan_tmpdir, ignore_errors=True)
        stan_tmpdir = tempfile.mkdtemp(prefix='stan_mv_mf_', dir=scratch)
        algorithm = 'meanfield'
        fit = model.variational(
            data=df_dict, seed=412, algorithm=algorithm,
            grad_samples=20, draws=1000,
            require_converged=False, show_console=True, refresh=1000,
            output_dir=stan_tmpdir
        )

    # Extract posterior draws via variational_sample_pd (matches original VIDRA).
    posteriors = fit.variational_sample_pd
    # Diagnostic: log draw count on first call to verify proper extraction
    if not hasattr(AS_multiVars, '_logged'):
        AS_multiVars._logged = True
        print(f"[DIAG] variational_sample_pd shape: {posteriors.shape}, "
              f"columns: {list(posteriors.columns[:8])}..., "
              f"slope unique values: {len(posteriors['slope'].unique())}/{ len(posteriors)}")
    slope_samples = posteriors['slope'].values
    intercept_samples = posteriors['intercept'].values
    sr_cols = [c for c in posteriors.columns if c.startswith('slope_random[')]
    slope_random_samples = posteriors[sorted(sr_cols)].values  # shape (draws, 5)

    # Build param_name -> samples dict for all slope/intercept parameters
    param_samples = {'slope': slope_samples, 'intercept': intercept_samples}
    for i in range(slope_random_samples.shape[1]):
        param_samples[f'slope_random[{i+1}]'] = slope_random_samples[:, i]

    # --- Meta-slope: pool slope_random posteriors for observed sources ---
    combination_observed = [
        list(pair) for pair in
        df[['GsourceLab', 'GqtlLab']].drop_duplicates().itertuples(
            index=False)
    ]
    slope_to_meta_analyse = []
    for combination in combination_observed:
        key = _get_slope_key(COMBINATION_SLOPE, combination)
        if key is not None:
            slope_to_meta_analyse.append(key)

    list_posteriors = []
    for slope_idx in slope_to_meta_analyse:
        col_name = f'slope_random[{slope_idx}]'
        if col_name in posteriors.columns:
            list_posteriors += posteriors[col_name].tolist()

    # --- Per-parameter posterior summaries ---
    rows = []
    for param_name, samples in param_samples.items():
        stats = clean_posteriorForAs(samples)
        stats['parameter'] = param_name
        rows.append(stats)

    # Add meta_slope row (pooled from observed slope_random components)
    if list_posteriors:
        meta_stats = clean_posteriorForAs(pd.Series(list_posteriors))
    else:
        # Fallback: use the fixed-effect slope if no source mapping matched
        meta_stats = clean_posteriorForAs(slope_samples)
    meta_stats['parameter'] = 'meta_slope'
    rows.append(meta_stats)

    # Build output DataFrame
    output = pd.DataFrame(rows)

    # Source/qtl: lists of unique values (matches original format)
    sources = str(sorted(set(int(x) for x in df['GsourceLab'])))
    qtls = str(sorted(set(int(x) for x in df['GqtlLab'])))

    output['gene'] = gene
    output['as_disease'] = str(df['as_disease'].iloc[0])
    output['n_variants'] = int(df['variant'].nunique())
    output['source'] = sources
    output['qtl'] = qtls
    output['model'] = f'multiple_variant_{algorithm}'

    # Clean up Stan temp dir
    shutil.rmtree(stan_tmpdir, ignore_errors=True)

    return output


def as_error_report(df, gene, error_msg=""):
    """Return a minimal error row when model fitting fails."""
    return pd.DataFrame([{
        'gene': gene,
        'as_disease': str(df['as_disease'].iloc[0]),
        'n_variants': int(df['variant'].nunique()),
        'model': 'error_report',
        'parameter': str(error_msg)[:500] if error_msg else 'error',
        'source': str(sorted(set(int(x) for x in df['GsourceLab']))),
        'qtl': str(sorted(set(int(x) for x in df['GqtlLab']))),
    }])


def fitASmodels(disease_df, gene, h1=0.1):
    """Decide single vs multi variant model and fit.

    Args:
        disease_df: DataFrame for ONE disease within a gene
        gene: gene identifier
        h1: hyperparameter

    Returns:
        DataFrame of posterior summaries
    """
    n_unique_variants = disease_df['variant'].nunique()
    try:
        if n_unique_variants > 1:
            result = AS_multiVars(disease_df, gene)
        else:
            result = AS_singleVars(disease_df, gene, h1)
        if result is None or result.empty:
            return as_error_report(disease_df, gene, 'empty result')
        return result
    except Exception as e:
        import traceback
        error_msg = f"{type(e).__name__}: {e} || {traceback.format_exc()[-400:]}"
        print(f"Stan error for {gene}/{disease_df['as_disease'].iloc[0]}: {error_msg}")
        return as_error_report(disease_df, gene, error_msg)


# =============================================================================
# Per-gene processing function (runs inside Spark UDF)
# =============================================================================

def process_gene(gene_df, h1=0.1):
    """Process all diseases for a single gene.

    Steps (matching original VIDRA_estimate_single_gene_linear_bayes.py):
      1. Fill missing annotation columns with defaults
      2. Dedup variants within (disease, source, qtl) — keep highest yc
      3. Filter out single-variant coding GWAS groups (GsourceLab==3)
      4. For each disease, run single or multi variant model
      5. Collect and return results

    Args:
        gene_df: pandas DataFrame for one gene (all diseases/variants)
        h1: hyperparameter for Stan models

    Returns:
        pandas DataFrame with posterior summaries per (disease, parameter)
    """
    if gene_df.empty:
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    gene = str(gene_df['as_gene'].iloc[0])

    # Verify models load before doing any processing
    try:
        get_models()
    except Exception as e:
        import traceback
        err = f"ModelLoad|{type(e).__name__}: {e}|{traceback.format_exc()[-500:]}"
        print(f"FATAL: Stan model load failed for {gene}: {err}")
        return as_error_report(gene_df.iloc[[0]], gene, err)

    # --- Step 1: Compute per-gene FoldX transform (as_foldx) ---
    # The FoldX ΔΔG score is normalised per-gene (relative to gene's max).
    # foldxDdq_raw holds the raw value from FoldX; we transform here.
    if 'foldxDdq_raw' in gene_df.columns:
        foldx_raw = pd.to_numeric(gene_df['foldxDdq_raw'], errors='coerce')
        max_score = foldx_raw.max()
        if pd.notna(max_score) and max_score != 0:
            inverted = max_score - foldx_raw
            max_inv = inverted.max()
            if max_inv > 0:
                gene_df['as_foldx'] = (1.0 - (inverted / max_inv)).fillna(1.0)
            else:
                gene_df['as_foldx'] = 1.0
        else:
            gene_df['as_foldx'] = 1.0

    # --- Fill annotation defaults (NON-missense-specific only) ---
    # Missense-specific columns (SIFT, PolyPhen, REVEL, AlphaMissense,
    # FoldX) arrive with NaN for non-missense variants — this is
    # intentional so per-gene mean imputation can work.  Only ensure
    # the columns exist; fill with NaN if absent.
    _missense_set = set(GENE_MEAN_IMPUTE_COLS)
    for col, default in ANNOTATION_DEFAULTS.items():
        if col in _missense_set:
            # Just ensure column exists (NaN placeholder); don't fill
            if col not in gene_df.columns:
                gene_df[col] = float('nan')
        else:
            _fill_annotation(gene_df, col, default)

    # --- Per-gene mean imputation for missense-specific annotations ---
    # For annotations only defined for missense variants (REVEL,
    # AlphaMissense, etc.), fill NaN with the per-gene mean of observed
    # values rather than a global benign default. This matches the
    # original pipeline (pre_processing_VIDRA_per_gene_pheno.py:440).
    for col in GENE_MEAN_IMPUTE_COLS:
        if col in gene_df.columns:
            gene_mean = gene_df[col].mean(skipna=True)
            if pd.notna(gene_mean):
                gene_df[col] = gene_df[col].fillna(gene_mean)

    # --- Hardcode LoF variants to maximally damaging (0.0) ---
    # stop_gained and start_lost variants should have all protein annotations
    # set to 0 (damaging) because missense-specific tools produce incorrect
    # benign defaults for them.
    # Matches original pipeline (pre_processing_VIDRA_per_gene_pheno.py:449).
    if 'most_severe_consequence' in gene_df.columns:
        lof_mask = gene_df['most_severe_consequence'].isin(LOF_CONSEQUENCES)
        if lof_mask.any():
            for col in GENE_MEAN_IMPUTE_COLS:
                if col in gene_df.columns:
                    gene_df.loc[lof_mask, col] = 0.0

    # --- Fill remaining NaN in gene-mean-imputed cols with global defaults ---
    # Catches genes where ALL variants had NaN (no observed values at all).
    for col in GENE_MEAN_IMPUTE_COLS:
        if col in gene_df.columns:
            gene_df[col] = gene_df[col].fillna(
                ANNOTATION_DEFAULTS.get(col, 1.0))

    # Replace any remaining NaN in numeric cols used by Stan
    numeric_cols = [
        'xc', 'xcse', 'yc', 'ycse', 'bO', 'bOse',
        'GsourceLab', 'GqtlLab'
    ]
    for col in numeric_cols:
        if col in gene_df.columns:
            gene_df[col] = gene_df[col].fillna(0.0)

    # Ensure integer types for source labels
    gene_df['GsourceLab'] = gene_df['GsourceLab'].astype(int)
    gene_df['GqtlLab'] = gene_df['GqtlLab'].astype(int)

    # --- Step 2: Dedup within (disease, source, qtl) groups ---
    # Sort by yc ascending → keep='last' retains highest effect size
    # (replicates filter_duplecate_lead without needing coloc_ref_dict.pkl)
    gene_df = gene_df.sort_values('yc')
    gene_df = gene_df.drop_duplicates(
        subset=['variant', 'as_disease', 'GsourceLab', 'GqtlLab'],
        keep='last'
    )

    # --- Step 3: Filter single-variant coding GWAS ---
    # Original: .filter(lambda x: not ((len(x['variant'])==1) & (x['GsourceLab']==3).all()))
    gene_df = gene_df.groupby(
        ['as_disease', 'GsourceLab', 'GqtlLab'], group_keys=False
    ).filter(
        lambda x: not (
            (len(x['variant']) == 1) and (x['GsourceLab'] == 3).all()
        )
    )

    if gene_df.empty:
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    # --- Step 4: Run model per disease ---
    results = gene_df.groupby('as_disease', group_keys=False).apply(
        lambda x: fitASmodels(x, gene, h1)
    )

    if results is None or results.empty:
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    # Ensure all output columns present (fill missing with None)
    for col in ALL_OUTPUT_COLS:
        if col not in results.columns:
            results[col] = None

    return results[ALL_OUTPUT_COLS]


# =============================================================================
# Spark output schema & UDF
# =============================================================================

result_schema = StructType([
    StructField("gene", StringType(), True),
    StructField("as_disease", StringType(), True),
    StructField("parameter", StringType(), True),
    StructField("model", StringType(), True),
    StructField("n_variants", IntegerType(), True),
    StructField("source", StringType(), True),
    StructField("qtl", StringType(), True),
    StructField("mean", DoubleType(), True),
    StructField("median", DoubleType(), True),
    StructField("pct_1", DoubleType(), True),
    StructField("pct_2_5", DoubleType(), True),
    StructField("pct_5", DoubleType(), True),
    StructField("pct_10", DoubleType(), True),
    StructField("pct_25", DoubleType(), True),
    StructField("pct_40", DoubleType(), True),
    StructField("pct_50", DoubleType(), True),
    StructField("pct_60", DoubleType(), True),
    StructField("pct_75", DoubleType(), True),
    StructField("pct_90", DoubleType(), True),
    StructField("pct_95", DoubleType(), True),
    StructField("pct_97_5", DoubleType(), True),
    StructField("pct_99", DoubleType(), True),
    StructField("pp_slope_pos", DoubleType(), True),
    StructField("pp_slope_neg", DoubleType(), True),
])


def main(args):
    spark = SparkSession.builder.appName("VIDRA_Bayesian_Analysis").getOrCreate()
    from pyspark.sql import functions as F

    ANALYSIS_PATH = f"gs://{args.bucket_name}/vidra_analysis_ready"
    ANNOTATIONS_PATH = f"gs://{args.bucket_name}/variant_annotations"
    OUTPUT_PATH = f"gs://{args.bucket_name}/vidra_results"
    h1 = float(args.h1)

    print("--- VIDRA Bayesian Analysis ---")
    print(f"Analysis data: {ANALYSIS_PATH}")
    print(f"Annotations:   {ANNOTATIONS_PATH}")
    print(f"Output:        {OUTPUT_PATH}")
    print(f"h1:            {h1}")

    # =========================================================================
    # 1. LOAD ANALYSIS-READY DATA (Hive-partitioned by as_gene, 1 file per gene)
    # =========================================================================
    # When --gene_list is provided, we use a literal .isin() filter which
    # guarantees Spark only lists+reads the requested gene partitions.
    # This is more reliable than a join-based filter for partition pruning.

    if args.gene_list:
        print(f"Loading target genes from {args.gene_list}")
        gene_list_df = spark.read.text(args.gene_list).withColumnRenamed("value", "as_gene")
        gene_list_df = gene_list_df.withColumn("as_gene", F.trim(F.col("as_gene")))
        # Collect to Python list for literal .isin() partition pruning
        target_genes = [row.as_gene for row in gene_list_df.collect()]
        target_genes = [g for g in target_genes if g]  # drop empty lines
        print(f"  {len(target_genes)} genes from gene list")
    else:
        target_genes = None

    # Read analysis-ready data with partition pruning
    df = spark.read.parquet(ANALYSIS_PATH)

    if target_genes:
        # Literal .isin() → Spark pushes this down to partition listing,
        # so only N directories are listed+read for N genes.
        df = df.filter(F.col("as_gene").isin(target_genes))

    # Test mode: add random genes to reach --test_genes total
    if args.test_mode:
        existing_genes = target_genes or []
        n_needed = max(0, args.test_genes - len(existing_genes))

        if n_needed > 0:
            print(f"TEST MODE: Adding {n_needed} random genes to reach {args.test_genes} total")
            all_genes_df = spark.read.parquet(ANALYSIS_PATH).select("as_gene").distinct()
            if existing_genes:
                pool = all_genes_df.filter(~F.col("as_gene").isin(existing_genes))
            else:
                pool = all_genes_df
            random_genes = [row.as_gene for row in
                            pool.orderBy(F.rand(seed=42)).limit(n_needed).collect()]
            all_target = existing_genes + random_genes
            print(f"  Total genes: {len(all_target)}")
            # Re-read with the expanded gene list for partition pruning
            df = spark.read.parquet(ANALYSIS_PATH).filter(
                F.col("as_gene").isin(all_target))

    print(f"Analysis data loaded: {df.count()} rows, "
          f"{df.select('as_gene').distinct().count()} genes")

    # =========================================================================
    # 2. JOIN VARIANT ANNOTATIONS (38 MB — Spark auto-broadcasts)
    # =========================================================================
    print(f"Loading variant annotations from {ANNOTATIONS_PATH} ...")
    annotations = spark.read.parquet(ANNOTATIONS_PATH)
    print(f"  Annotations: {annotations.count()} rows")

    # Drop any annotation columns from base data that will come from the join.
    # Step 1 output now only carries as_clinicalSignificance (from OT ClinVar
    # evidence); the guard handles backwards compatibility if older Step 1
    # output still contains annotation columns.
    cols_to_drop = [c for c in ANNOTATION_COLS
                    if c in df.columns and c != 'as_clinicalSignificance']
    if 'as_clinicalSignificance' in df.columns:
        df = df.withColumnRenamed('as_clinicalSignificance',
                                 '_base_as_clinicalSignificance')
    if cols_to_drop:
        df = df.drop(*cols_to_drop)

    # Left-join on variant
    df = df.join(annotations, on="variant", how="left")
    print(f"  After join: {df.count()} rows")

    # Merge as_clinicalSignificance: prefer OT ClinVar over VEP
    if '_base_as_clinicalSignificance' in df.columns:
        df = df.withColumn(
            'as_clinicalSignificance',
            F.when(
                F.coalesce(F.col('_base_as_clinicalSignificance'), F.lit(0.0)) > 0.0,
                F.col('_base_as_clinicalSignificance')
            ).otherwise(
                F.coalesce(F.col('as_clinicalSignificance'), F.lit(0.0))
            )
        ).drop('_base_as_clinicalSignificance')

    # =========================================================================
    # 3. APPLY DIRECTION INVERSIONS & NORMALISATIONS
    # =========================================================================
    # protein_prior convention: 0 = damaging, 1 = benign.
    # These annotations have raw high = more damaging, so invert: 1 - x
    for col in INVERSION_COLS:
        if col in df.columns:
            df = df.withColumn(col, F.lit(1.0) - F.col(col))
    print(f"Applied 1-x inversions to: {', '.join(c for c in INVERSION_COLS if c in df.columns)}")

    # Normalise & invert as_consequence: ordinal 0-12 → 1 - code/12
    # Convention: 0 = most damaging (high consequence code), 1 = benign (low code)
    # Currently commented out in both Stan models but kept correct for future use.
    _N_CONSEQUENCE_CATS = 12  # max ordinal code in CONSEQUENCE_CATEGORIES
    if 'as_consequence' in df.columns:
        df = df.withColumn('as_consequence',
                          F.lit(1.0) - F.col('as_consequence').cast('double') / _N_CONSEQUENCE_CATS)
        print(f"Normalised & inverted as_consequence to [0, 1] (max code={_N_CONSEQUENCE_CATS})")

    # Normalise as_clinicalSignificance to [0, 1]
    if 'as_clinicalSignificance' in df.columns:
        cs_max = df.agg(F.max('as_clinicalSignificance')).collect()[0][0]
        if cs_max is not None and cs_max > 1.0:
            print(f"Rescaling as_clinicalSignificance (max={cs_max:.2f}) to [0, 1]")
            df = df.withColumn('as_clinicalSignificance',
                              F.col('as_clinicalSignificance') / float(cs_max))
        else:
            print(f"as_clinicalSignificance already in [0, 1] (max={cs_max})")

    # Fill NAs with defaults — but EXCLUDE missense-specific columns.
    # Those must remain NaN so that per-gene mean imputation in
    # process_gene() can distinguish observed values from missing ones
    # (matching original pre_processing_VIDRA_per_gene_pheno.py logic).
    _skip_spark_fill = set(GENE_MEAN_IMPUTE_COLS)
    fill_dict = {k: v for k, v in ANNOTATION_DEFAULTS.items()
                 if k in df.columns and isinstance(v, (int, float))
                 and not (isinstance(v, float) and v != v)  # skip NaN
                 and k not in _skip_spark_fill}
    df = df.fillna(fill_dict)
    print(f"Filled NA defaults for {len(fill_dict)} annotation columns "
          f"(deferred {len(_skip_spark_fill)} missense-specific cols for per-gene imputation)")

    # =========================================================================
    # 4. REPARTITION & RUN STAN MODELS
    # =========================================================================
    n_genes = df.select("as_gene").distinct().count()
    total_rows = df.count()
    print(f"Ready for Stan: {total_rows} rows across {n_genes} genes")

    # Repartition by gene so each gene gets its own Spark task
    df = df.repartition(n_genes, "as_gene")

    # Run per-gene analysis via applyInPandas (Spark 3.x)
    results = df.groupby("as_gene").applyInPandas(
        lambda pdf: process_gene(pdf, h1=h1),
        schema=result_schema
    )

    results.write.mode("overwrite").parquet(OUTPUT_PATH)

    n_result_rows = spark.read.parquet(OUTPUT_PATH).count()
    print(f"Analysis complete. {n_result_rows} result rows at {OUTPUT_PATH}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="VIDRA Bayesian Analysis on Dataproc Serverless")
    parser.add_argument("--bucket_name", required=True,
                        help="GCS bucket name (without gs://)")
    parser.add_argument("--h1", type=float, default=0.1,
                        help="Hyperparameter h1 (default: 0.1)")
    parser.add_argument("--test_mode", action="store_true",
                        help="Run on a random subset of genes")
    parser.add_argument("--test_genes", type=int, default=50,
                        help="Number of genes in test mode (default: 50)")
    parser.add_argument("--gene_list", type=str,
                        help="Path to a text file (e.g., gs://bucket/genes.txt) containing ENSG IDs (one per line) to filter the analysis")
    args = parser.parse_args()
    main(args)
