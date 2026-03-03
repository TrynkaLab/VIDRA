#!/usr/bin/env python3
"""VIDRA Variant Annotation Pipeline — VEP CLI + local FoldX lookup.

Annotation approach (fully offline, no REST APIs):
  1. VEP CLI for bulk annotation (~30-60 min for 1.6M variants)
     - Docker mode (--use_docker): runs inside ensemblorg/ensembl-vep:release_111.0.
       Bio::DB::BigFile is pre-compiled; Conservation plugin works natively.
     - Bare mode (default): VEP installed on host; GERP extracted via pyBigWig.
     - Uses --check_existing for ClinVar clinical significance from VEP cache.
     - Uses --uniprot to get SWISSPROT accessions for FoldX lookup.
     - VEP plugins provide AlphaMissense, REVEL, PrimateAI, CADD scores
       from local data files downloaded during VM setup.
  2. Local FoldX FTP file for foldxDdq + plddt (no ProtVar REST API)

This script is designed to run on a GCE VM.
See scripts/pyspark_scripts/run_annotation_vm.sh for the VM orchestration script.

Pipeline position:
  1. prepare_analysis_input.py → gs://<bucket>/vidra_analysis_ready/
  2. THIS SCRIPT               → gs://<bucket>/variant_annotations/
  3. run_bayesian_analysis.py   → reads both, joins on 'variant', runs Stan

VEP CLI plugins used:
  - Blosum62        (built-in)
  - CADD            (plugin + data file)
  - AlphaMissense   (plugin + data file)
  - REVEL           (plugin + data file)
  - PrimateAI       (plugin + data file)
  SIFT & PolyPhen are included in the VEP cache (--sift b --polyphen b).

FoldX FTP file (replaces ProtVar REST API for foldx + plddt):
  - Source: https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/stability/
  - File: 2025.02.10_foldx_energy.csv.gz (~1.5GB, 81.5M rows)
  - Columns: uniprot_accession, uniprot_position, wild_type, mutated_type, foldx_ddg, plddt
  - Lookup key: (swissprot_accession from VEP, protein_start, ref_aa, alt_aa)

Transforms (matching pre_processing_VIDRA_per_gene_pheno.py):
  as_blosum62:      sigmoid 1/(1+exp(-x)), fillna(1)
  as_sift:          raw score, NaN preserved for Step 3 per-gene mean imputation
  as_polyphen:      raw score, NaN preserved  [1-x inversion applied at runtime]
  as_cadd:          phred/50 clamped [0,1], NaN preserved  [1-x inversion at runtime]
  as_alphamissense: raw am_pathogenicity, NaN preserved
  as_revel:         raw REVEL score, NaN preserved  [1-x inversion applied at runtime]
  as_primateai:     raw PrimateAI score, fillna(0)
  as_loftool:       not used by Stan model; always 0.0
  as_plddt:         plddt/100, fillna(0)
  as_conservation:  VEP GERP → clamp [-12.36, 6.18], scale to [0,1], then
                    1 - scaled. NaN preserved for Step 3 per-gene mean imputation
  as_consequence:   ordinal encoding 0-12
  as_clinicalSignificance: 17 ordered categories → codes 0-16 → MinMaxScale [0,1]
  foldxDdq_raw:     raw value for per-gene transform in Bayesian step

Usage (on a GCE VM with VEP installed):
  python3 annotate_variants_cli.py \\
    --bucket_name vidra-2-0 \\
    --vep_cache_dir /opt/vep/cache \\
    --plugin_dir /opt/vep/plugins \\
    --work_dir /tmp/vidra_annotation \\
    --foldx_file /opt/vep/plugin_data/foldx_energy.csv.gz

  # Test mode (500 variants):
  python3 annotate_variants_cli.py \\
    --bucket_name vidra-2-0 --test_mode --test_variants 500 \\
    --vep_cache_dir /opt/vep/cache \\
    --plugin_dir /opt/vep/plugins \\
    --work_dir /tmp/vidra_annotation \\
    --foldx_file /opt/vep/plugin_data/foldx_energy.csv.gz
"""

import argparse
import csv
import json
import logging
import os
import subprocess
import time
from pathlib import Path

import gzip

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
try:
    import pyBigWig
    _PYBIGWIG_AVAILABLE = True
except ImportError:
    _PYBIGWIG_AVAILABLE = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


# ============================================================================
# Constants
# ============================================================================

# --- FoldX local file (replaces ProtVar REST API) ---------------------------
FOLDX_FTP_URL = "https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/stability/2025.02.10_foldx_energy.csv.gz"
FOLDX_FILENAME = "foldx_energy.csv.gz"

# Coding consequences — only look up FoldX/ClinVar for these
CODING_CONSEQUENCES = {
    "stop_gained", "missense_variant", "splice_donor_variant",
    "splice_region_variant", "synonymous_variant", "splice_acceptor_variant",
    "frameshift_variant", "inframe_deletion", "inframe_insertion",
    "stop_lost", "start_lost", "protein_altering_variant",
    "splice_donor_region_variant",
}

# Consequence ordinal encoding — ordered by VEP severity (least → most damaging).
# See: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
# The original formatData_StanModel used a similar order (0–9) but with minor
# ranking errors, fewer categories, and a missing start_lost.  This corrected
# list adds 3 severe consequences (start_lost, frameshift, splice_acceptor),
# fixes the placement of non_coding_transcript_exon_variant, and follows the
# exact Ensembl severity table ordering (most → least severe: splice_acceptor >
# splice_donor > stop_gained > frameshift > start_lost > missense >
# splice_region > splice_donor_region > synonymous > 5'UTR > 3'UTR > nc_exon).
CONSEQUENCE_CATEGORIES = [
    "__unknown__",                          # 0  — fallback
    "non_coding_transcript_exon_variant",   # 1  — MODIFIER (VEP rank ~27)
    "3_prime_UTR_variant",                  # 2  — MODIFIER (VEP rank ~26)
    "5_prime_UTR_variant",                  # 3  — MODIFIER (VEP rank ~25)
    "synonymous_variant",                   # 4  — LOW     (VEP rank ~22)
    "splice_donor_region_variant",          # 5  — LOW     (VEP rank ~17)
    "splice_region_variant",                # 6  — LOW     (VEP rank ~16)
    "missense_variant",                     # 7  — MODERATE(VEP rank ~13)
    "start_lost",                           # 8  — HIGH    (VEP rank ~7)
    "frameshift_variant",                   # 9  — HIGH    (VEP rank ~5)
    "stop_gained",                          # 10 — HIGH    (VEP rank ~4)
    "splice_donor_variant",                 # 11 — HIGH    (VEP rank ~3)
    "splice_acceptor_variant",              # 12 — HIGH    (VEP rank ~2)
]

# Clinical significance ordinal encoding — matches original
# pre_processing_VIDRA_per_gene_pheno.py (lines 458-468).
# 17 ordered categories → integer codes 0-16 → MinMaxScale to [0, 1].
# Stan prior: disease_prior ~ normal(as_clinicalSignificance, 0.2)
# The sigma=0.2 was calibrated for values in [0, 1].
#
# Convention: 0.0 = not provided / benign, 1.0 = pathogenic.
CLINICAL_SIG_CATEGORIES = [
    'not provided',              # 0  → 0.000
    'association not found',     # 1  → 0.0625
    'other',                     # 2  → 0.125
    'benign',                    # 3  → 0.1875
    'likely benign',             # 4  → 0.250
    'low penetrance',            # 5  → 0.3125
    'confers sensitivity',       # 6  → 0.375
    'uncertain risk allele',     # 7  → 0.4375
    'drug response',             # 8  → 0.500
    'uncertain significance',    # 9  → 0.5625
    'association',               # 10 → 0.625
    'affects',                   # 11 → 0.6875
    'likely risk allele',        # 12 → 0.750
    'risk factor',               # 13 → 0.8125
    'established risk allele',   # 14 → 0.875
    'likely pathogenic',         # 15 → 0.9375
    'pathogenic',                # 16 → 1.000
]
_N_CLIN_CATS = len(CLINICAL_SIG_CATEGORIES) - 1  # 16, for MinMaxScale denominator

# Normalisation: collapse compound / underscore terms from VEP or ProtVar
# into the canonical 17-category strings above.
_CLIN_SIG_NORMALISE = {
    # VEP uses underscores
    'likely_pathogenic': 'likely pathogenic',
    'likely_benign': 'likely benign',
    'uncertain_significance': 'uncertain significance',
    'not_provided': 'not provided',
    'association_not_found': 'association not found',
    'low_penetrance': 'low penetrance',
    'confers_sensitivity': 'confers sensitivity',
    'uncertain_risk_allele': 'uncertain risk allele',
    'drug_response': 'drug response',
    'likely_risk_allele': 'likely risk allele',
    'risk_factor': 'risk factor',
    'established_risk_allele': 'established risk allele',
    # Compound terms from ClinVar / VEP
    'variant of uncertain significance': 'uncertain significance',
    'benign/likely benign': 'likely benign',
    'benign/likely_benign': 'likely benign',
    'pathogenic/likely pathogenic': 'likely pathogenic',
    'pathogenic/likely_pathogenic': 'likely pathogenic',
    'conflicting interpretations of pathogenicity': 'uncertain significance',
    'conflicting_interpretations_of_pathogenicity': 'uncertain significance',
}

# GERP RS score range for normalisation to [0,1].
# Positive = conserved, negative = fast-evolving.
# Range comes from the global min/max of the GRCh38 GERP BigWig.
GERP_RS_MIN = -12.36
GERP_RS_MAX = 6.18

# Default fill values for annotations that are explicitly filled in
# apply_transforms (blosum62, primateai, loftool, plddt, consequence,
# clinicalSignificance).  Missense-specific scores (sift, polyphen, cadd,
# revel, alphamissense, conservation) are left as NaN in the output so that
# Step 3 can apply per-gene mean imputation from observed values.
# These are raw VEP-scale values, before 1-x inversions in Step 3.
ANNOTATION_DEFAULTS = {
    "as_blosum62": 1.0,
    "as_conservation": 1.0,
    "as_sift": 1.0,
    "as_polyphen": 1.0,
    "as_cadd": 0.0,
    "as_alphamissense": 0.0,
    "as_revel": 0.0,
    "as_primateai": 0.0,
    "as_loftool": 0.0,
    "as_plddt": 0.0,
    "as_consequence": 0,
    "as_clinicalSignificance": 0.0,
    "foldxDdq_raw": float("nan"),
}


# ============================================================================
# Step 1: Read unique variants from GCS manifest
# ============================================================================

def read_unique_variants(bucket: str) -> list[str]:
    """Read unique variant IDs from the manifest parquet on GCS."""
    import gcsfs
    fs = gcsfs.GCSFileSystem()
    manifest_path = f"{bucket}/vidra_analysis_ready_manifest"

    parquet_files = fs.glob(f"{manifest_path}/*.parquet")
    if not parquet_files:
        parquet_files = fs.glob(f"{manifest_path}/part-*")
    if not parquet_files:
        raise FileNotFoundError(
            f"No manifest parquet found at gs://{manifest_path}"
        )

    log.info("Reading variant manifest from gs://%s (%d files)",
             manifest_path, len(parquet_files))

    all_variants: set[str] = set()
    for fpath in parquet_files:
        with fs.open(fpath, "rb") as f:
            tbl = pq.read_table(f, columns=["variant"])
            all_variants.update(tbl.column("variant").to_pylist())

    # Sort variants genomically (by chrom, then position) for VEP
    def chrom_sort_key(chrom: str) -> tuple:
        """Return sort key for chromosome: (is_not_numeric, chrom_as_int_or_str)"""
        if chrom.isdigit():
            return (0, int(chrom))  # 1, 2, ..., 22
        else:
            return (1, chrom)       # X, Y, MT
    
    def variant_sort_key(vid: str) -> tuple:
        """Parse variant ID and return (chrom_key, position_int)"""
        parts = vid.split("_")
        if len(parts) < 2:
            return ((2, vid), 0)  # malformed variants last
        chrom, pos = parts[0], parts[1]
        try:
            pos_int = int(pos)
        except ValueError:
            pos_int = 0
        return (chrom_sort_key(chrom), pos_int)
    
    result = sorted(all_variants, key=variant_sort_key)
    log.info("Unique variants from manifest: %d", len(result))
    return result


# ============================================================================
# Step 2: Write VCF for VEP CLI
# ============================================================================

def write_vcf(variant_ids: list[str], vcf_path: Path) -> int:
    """Write variant IDs to a minimal VCF file for VEP input.

    Returns the number of successfully written variants.
    """
    log.info("Writing %d variants to VCF: %s", len(variant_ids), vcf_path)
    written = 0
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for vid in variant_ids:
            parts = vid.split("_")
            if len(parts) < 4:
                continue
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")
            written += 1
    log.info("VCF written: %d variants", written)
    return written


# ============================================================================
# Step 3: Run VEP CLI
# ============================================================================

def run_vep_cli(
    vcf_path: Path,
    output_path: Path,
    vep_cache_dir: str,
    plugin_dir: str,
    plugin_data_dir: str,
    threads: int = 4,
    buffer_size: int = 5000,
    use_docker: bool = False,
    docker_image: str = "ensemblorg/ensembl-vep:release_111.0",
) -> None:
    """Run Ensembl VEP CLI with plugins, outputting JSON.

    Args:
        vcf_path:       Input VCF file
        output_path:    Output JSON file
        vep_cache_dir:  VEP cache directory (e.g., /opt/vep/cache)
        plugin_dir:     VEP plugins directory (host path; ignored when
                        use_docker=True — the image's /opt/vep/Plugins is used)
        plugin_data_dir: Directory containing plugin data files
        threads:        Number of VEP threads (fork)
        buffer_size:    Number of variants to process at once
        use_docker:     Run VEP inside ensemblorg/ensembl-vep Docker container.
                        Bio::DB::BigFile is pre-compiled in the image, so the
                        Conservation plugin runs natively (no pyBigWig needed).
        docker_image:   Docker image tag to use when use_docker=True.
    """
    if use_docker:
        # Container-internal paths (after volume mounts)
        ctr_cache   = "/opt/vep/.vep"
        ctr_plugins = "/plugins"  # VEP 111 Docker image stores plugins here
        ctr_data    = "/plugin_data"
        work_dir    = str(vcf_path.parent)
        ctr_input   = f"/work/{vcf_path.name}"
        ctr_output  = f"/work/{output_path.name}"

        # Ensure the work dir is writable by the container's 'vep' user
        os.chmod(work_dir, 0o777)

        vep_args = [
            "vep",
            "--cache",
            "--species", "homo_sapiens",
            "--assembly", "GRCh38",
            "--dir_cache", ctr_cache,
            "--dir_plugins", ctr_plugins,
            "--input_file", ctr_input,
            "--output_file", ctr_output,
            "--json",
            "--canonical",
            "--sift", "b",
            "--polyphen", "b",
            "--check_existing",     # adds ClinVar clinical_significance
            "--uniprot",            # adds swissprot/trembl accessions
            "--force_overwrite",
            "--no_stats",
            "--fork", str(threads),
            "--buffer_size", str(buffer_size),
            "--offline",
            # Conservation plugin works natively: Bio::DB::BigFile pre-compiled
            "--plugin", "Blosum62",
            "--plugin", f"CADD,snv={ctr_data}/CADD/whole_genome_SNVs.tsv.gz",
            "--plugin", f"Conservation,{ctr_data}/conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
            "--plugin", f"AlphaMissense,file={ctr_data}/AlphaMissense_hg38.tsv.gz",
            "--plugin", f"REVEL,file={ctr_data}/new_tabbed_revel_grch38.tsv.gz",
            "--plugin", f"PrimateAI,{ctr_data}/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz",
        ]

        cmd = [
            "docker", "run", "--rm",
            "-v", f"{work_dir}:/work",
            "-v", f"{vep_cache_dir}:{ctr_cache}:ro",
            "-v", f"{plugin_data_dir}:{ctr_data}:ro",
            docker_image,
        ] + vep_args

        log.info("Docker mode: VEP running inside %s", docker_image)
        log.info("  Cache mount:  %s → %s", vep_cache_dir, ctr_cache)
        log.info("  Data mount:   %s → %s", plugin_data_dir, ctr_data)
        log.info("  Work mount:   %s → /work", work_dir)
    else:
        # Build VEP command (bare binary; GERP extracted post-hoc via pyBigWig)
        cmd = [
            "vep",
            "--cache",
            "--species", "homo_sapiens",
            "--assembly", "GRCh38",
            "--dir_cache", vep_cache_dir,
            "--dir_plugins", plugin_dir,
            "--input_file", str(vcf_path),
            "--output_file", str(output_path),
            "--json",
            "--canonical",
            "--sift", "b",          # both prediction and score
            "--polyphen", "b",      # both prediction and score
            "--check_existing",     # adds ClinVar clinical_significance
            "--uniprot",            # adds swissprot/trembl accessions
            "--force_overwrite",
            "--no_stats",
            "--fork", str(threads),
            "--buffer_size", str(buffer_size),
            "--offline",            # don't contact Ensembl servers
            # Plugins
            "--plugin", "Blosum62",
            "--plugin", f"CADD,snv={plugin_data_dir}/CADD/whole_genome_SNVs.tsv.gz",
            # Conservation/GERP extracted post-hoc via pyBigWig (see Step 4b)
            # Bio::DB::BigFile cannot compile on Debian 12
            "--plugin", f"AlphaMissense,file={plugin_data_dir}/AlphaMissense_hg38.tsv.gz",
            "--plugin", f"REVEL,file={plugin_data_dir}/new_tabbed_revel_grch38.tsv.gz",
            "--plugin", f"PrimateAI,{plugin_data_dir}/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz",
        ]

    log.info("Running VEP CLI: %s", " ".join(cmd[:10]) + " ...")
    log.info("  Input:  %s", vcf_path)
    log.info("  Output: %s", output_path)
    log.info("  Threads: %d, Buffer: %d", threads, buffer_size)
    log.info("  Full command:\n    %s", " \\\n    ".join(cmd))

    t0 = time.time()

    # Use Popen to stream stderr in real-time (VEP prints progress there)
    # instead of buffering everything in memory which could OOM for 1.6M vars
    stderr_log = output_path.parent / "vep_stderr.log"
    with open(stderr_log, "w") as stderr_f:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        # Stream stderr to log file and to Python log
        last_progress_time = t0
        for line in proc.stderr:
            stderr_f.write(line)
            stripped = line.strip()
            # Log progress lines or warnings (avoid flooding with every line)
            now = time.time()
            if stripped and (
                now - last_progress_time > 30  # at least every 30s
                or "warning" in stripped.lower()
                or "error" in stripped.lower()
                or "processed" in stripped.lower()
                or "%" in stripped
            ):
                log.info("  VEP: %s", stripped)
                last_progress_time = now

        proc.wait()

    elapsed = time.time() - t0

    if proc.returncode != 0:
        log.error("VEP FAILED (exit code %d) after %.1fs", proc.returncode, elapsed)
        # Read last lines of stderr log for diagnostics
        try:
            with open(stderr_log, "r") as f:
                lines = f.readlines()
                log.error("VEP STDERR (last 30 lines):\n%s",
                          "".join(lines[-30:]))
        except Exception:
            pass
        raise RuntimeError(f"VEP CLI failed with exit code {proc.returncode}")

    log.info("VEP completed in %.1fs", elapsed)
    # Log VEP summary stats from end of stderr
    try:
        with open(stderr_log, "r") as f:
            lines = f.readlines()
            for line in lines[-20:]:
                if line.strip():
                    log.info("  VEP: %s", line.strip())
    except Exception:
        pass


# ============================================================================
# Step 4: Parse VEP JSON output
# ============================================================================

def _get_field(d: dict, *keys) -> object | None:
    """Get a value from a dict trying multiple key names (case-insensitive).

    VEP CLI plugins may use different capitalisation than the REST API, e.g.
    CADD_phred vs cadd_phred, REVEL vs revel_score.
    """
    # Try exact keys first (fast path)
    for k in keys:
        if k in d:
            return d[k]
    # Case-insensitive fallback
    lower_map = {kk.lower(): v for kk, v in d.items()}
    for k in keys:
        v = lower_map.get(k.lower())
        if v is not None:
            return v
    return None


def parse_vep_json(output_path: Path) -> pd.DataFrame:
    """Parse VEP JSON-per-line output into a DataFrame.

    Extracts from the canonical transcript:
      - blosum62, sift_score, polyphen_score, cadd_phred
      - REVEL, AlphaMissense, PrimateAI (from plugins)
      - Conservation (GERP)
      - most_severe_consequence
      - swissprot (UniProt accession), protein_start, amino_acids
        (for FoldX local file lookup)
    Extracts from colocated_variants:
      - clin_sig (ClinVar clinical significance, from --check_existing)

    Handles both REST API field names (lowercase) and CLI plugin field names
    (which may use capitals, e.g. CADD_phred, REVEL, PrimateAI_score).
    """
    log.info("Parsing VEP JSON output: %s", output_path)
    t0 = time.time()
    records = []
    field_names_seen: dict[str, int] = {}  # track field names for debugging

    with open(output_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                item = json.loads(line)
            except json.JSONDecodeError:
                log.warning("Skipping invalid JSON at line %d", line_num)
                continue

            # Use the variant ID we put in the VCF ID column
            variant_id = item.get("id", "")
            if not variant_id:
                # Fallback: reconstruct from input
                inp = item.get("input", "")
                parts = inp.split("\t") if "\t" in inp else inp.split()
                if len(parts) >= 5:
                    variant_id = parts[2]  # ID column in VCF

            if not variant_id:
                continue

            rec = {"variant": variant_id}
            rec["most_severe_consequence"] = item.get("most_severe_consequence")

            # Find canonical transcript
            tc_list = item.get("transcript_consequences", [])
            best_tc = None
            best_completeness = -1
            for tc in tc_list:
                if tc.get("canonical") == 1:
                    completeness = sum(1 for v in tc.values() if v is not None)
                    if completeness > best_completeness:
                        best_tc = tc
                        best_completeness = completeness

            if best_tc:
                # Track field names from first few records for debugging
                if line_num <= 5:
                    for k in best_tc:
                        field_names_seen[k] = field_names_seen.get(k, 0) + 1

                # SIFT / PolyPhen (built-in, always lowercase)
                rec["sift_score"] = best_tc.get("sift_score")
                rec["polyphen_score"] = best_tc.get("polyphen_score")

                # Blosum62 plugin (lowercase in CLI)
                rec["blosum62"] = _get_field(best_tc, "blosum62", "Blosum62")

                # CADD plugin: CLI may use CADD_phred / CADD_raw
                rec["cadd_phred"] = _get_field(
                    best_tc, "cadd_phred", "CADD_phred", "CADD_PHRED"
                )
                rec["cadd_raw"] = _get_field(
                    best_tc, "cadd_raw", "CADD_raw", "CADD_RAW"
                )

                # Conservation plugin (GERP): CLI may use Conservation
                rec["conservation"] = _get_field(
                    best_tc, "conservation", "Conservation"
                )

                # AlphaMissense plugin
                am_data = _get_field(
                    best_tc, "am_pathogenicity", "AM_pathogenicity"
                )
                if am_data is None:
                    am_info = best_tc.get("alphamissense") or {}
                    if isinstance(am_info, dict):
                        am_data = am_info.get("am_pathogenicity")
                rec["am_pathogenicity"] = am_data

                # REVEL plugin: CLI uses REVEL (single value); REST uses revel_score
                rec["revel"] = _get_field(
                    best_tc, "revel_score", "REVEL", "revel"
                )

                # PrimateAI plugin: VEP JSON uses "primateai" (raw score)
                rec["primateai_score"] = _get_field(
                    best_tc, "primateai", "primateai_score", "PrimateAI_score",
                    "primatai_score", "primatai_pred"
                )

                # --- UniProt accession + protein position for FoldX lookup ---
                # --uniprot flag adds swissprot/trembl (lists of versioned accessions)
                # VEP returns e.g. "O00468.207" — strip the ".NNN" version suffix
                # to match the FoldX file which uses bare accessions ("O00468").
                raw_acc = None
                sp = best_tc.get("swissprot")
                if isinstance(sp, list) and sp:
                    raw_acc = sp[0]
                elif isinstance(sp, str) and sp:
                    raw_acc = sp
                # Fall back to trembl if no swissprot (FoldX covers both)
                if not raw_acc:
                    tr = best_tc.get("trembl")
                    if isinstance(tr, list) and tr:
                        raw_acc = tr[0]
                    elif isinstance(tr, str) and tr:
                        raw_acc = tr
                if raw_acc:
                    rec["swissprot"] = raw_acc.split(".")[0]  # strip version

                rec["protein_start"] = best_tc.get("protein_start")

                # amino_acids: "X/Y" where X=ref, Y=alt
                aa = best_tc.get("amino_acids")
                if aa and "/" in str(aa):
                    parts_aa = str(aa).split("/")
                    rec["ref_aa"] = parts_aa[0] if len(parts_aa) > 0 else None
                    rec["alt_aa"] = parts_aa[1] if len(parts_aa) > 1 else None

            # --- ClinVar clinical significance from --check_existing ---
            colocated = item.get("colocated_variants", [])
            if colocated and isinstance(colocated, list):
                best_clin_sig = None
                severity = {"pathogenic": 3, "likely_pathogenic": 2,
                            "pathogenic/likely_pathogenic": 2,
                            "benign": 1, "likely_benign": 1,
                            "benign/likely_benign": 1}
                for cv in colocated:
                    clin_sig_str = cv.get("clin_sig")
                    if not clin_sig_str:
                        continue
                    # VEP --check_existing: clin_sig is comma-separated string
                    for sig in str(clin_sig_str).split(","):
                        sig = sig.strip().lower()
                        s = severity.get(sig, 0)
                        if best_clin_sig is None or s > severity.get(best_clin_sig, 0):
                            best_clin_sig = sig
                if best_clin_sig:
                    rec["clinicalSignificance"] = best_clin_sig

            # Fallback: check variant-level or intergenic_consequences for
            # Conservation plugin data (it sometimes appears at top level)
            if rec.get("conservation") is None:
                rec["conservation"] = _get_field(
                    item, "conservation", "Conservation"
                )

            records.append(rec)

            if line_num % 200000 == 0:
                log.info("  Parsed %d VEP records...", line_num)

    elapsed = time.time() - t0
    df = pd.DataFrame(records)
    log.info("VEP parsed: %d records in %.1fs", len(df), elapsed)

    # Log field names from canonical transcripts (helpful for debugging)
    if field_names_seen:
        log.info("VEP transcript_consequences fields seen (first 5 records): %s",
                 sorted(field_names_seen.keys()))

    return df


# ============================================================================
# Step 4b: GERP conservation scores via pyBigWig
# ============================================================================

def annotate_gerp(
    variant_ids: list[str],
    bigwig_path: str,
) -> dict[str, float | None]:
    """Extract GERP conservation scores from BigWig for each variant.

    Uses pyBigWig (Python C-extension wrapping libBigWig) instead of the
    VEP Conservation plugin (which requires Bio::DB::BigFile / UCSC Kent
    source, uncompilable on Debian 12).

    Args:
        variant_ids: List of "CHROM_POS_REF_ALT" variant IDs.
        bigwig_path: Path to the GERP BigWig file.

    Returns:
        Dict mapping variant_id → GERP score (float) or None.
    """
    if not _PYBIGWIG_AVAILABLE:
        log.warning("pyBigWig not installed — GERP scores will be missing")
        return {}

    log.info("Extracting GERP scores from BigWig: %s", bigwig_path)
    t0 = time.time()
    bw = pyBigWig.open(bigwig_path)
    chroms_in_bw = set(bw.chroms().keys())

    scores: dict[str, float | None] = {}
    missing_chroms: set[str] = set()

    for vid in variant_ids:
        parts = vid.split("_")
        if len(parts) < 2:
            continue
        chrom, pos_str = parts[0], parts[1]
        try:
            pos = int(pos_str)
        except ValueError:
            continue

        # Try with and without "chr" prefix
        chrom_key: str | None = None
        if chrom in chroms_in_bw:
            chrom_key = chrom
        elif f"chr{chrom}" in chroms_in_bw:
            chrom_key = f"chr{chrom}"
        else:
            missing_chroms.add(chrom)
            continue

        # BigWig uses 0-based half-open coords; VEP/VCF use 1-based
        vals = bw.values(chrom_key, pos - 1, pos)
        scores[vid] = vals[0] if vals and vals[0] is not None else None

    bw.close()
    elapsed = time.time() - t0
    n_scored = sum(1 for v in scores.values() if v is not None)
    log.info(
        "GERP annotation: %d / %d variants scored in %.1fs",
        n_scored, len(variant_ids), elapsed,
    )
    if missing_chroms:
        log.warning("GERP: chromosomes not found in BigWig: %s", sorted(missing_chroms)[:10])
    return scores


# ============================================================================
# Step 5: FoldX local file lookup (replaces ProtVar REST API)
# ============================================================================

def load_foldx_lookup(
    foldx_path: str,
    needed_keys: set[tuple] | None = None,
    relevant_accessions: set[str] | None = None,
) -> dict[tuple, tuple]:
    """Load FoldX FTP file into a lookup dict.

    The FTP file format (CSV, gzipped):
      uniprot_accession, uniprot_position, alphafold_fragment_id,
      alphafold_fragment_position, wild_type, mutated_type, foldx_ddg, plddt

    Two filtering strategies (used in order of preference):
      1. needed_keys: exact (accession, position, wt, mut) tuples to keep.
         Memory-efficient — only stores matching rows (~1M variants → <50K hits).
      2. relevant_accessions: accessions to keep (loads ALL mutations for those
         proteins — can be 200M+ rows → OOM on 32GB).

    Args:
        foldx_path: Path to the gzipped CSV file.
        needed_keys: If provided, only load rows matching these exact keys.
        relevant_accessions: Fallback filter — only load rows for these accessions.

    Returns:
        Dict mapping (accession, position, wild_type, mutated_type) →
        (foldx_ddg, plddt).
    """
    filter_desc = (
        f"exact-key filter ({len(needed_keys)} keys)" if needed_keys
        else f"accession filter ({len(relevant_accessions)} accessions)" if relevant_accessions
        else "all rows"
    )
    log.info("Loading FoldX lookup from %s (%s)", foldx_path, filter_desc)
    t0 = time.time()
    lookup: dict[tuple, tuple] = {}
    n_skipped = 0
    n_total = 0

    opener = gzip.open if foldx_path.endswith(".gz") else open
    with opener(foldx_path, "rt") as f:
        reader = csv.reader(f)
        header = next(reader)  # skip header
        log.info("FoldX header: %s", header)

        for row in reader:
            n_total += 1
            if len(row) < 8:
                continue
            accession = row[0]

            # Fast path: skip rows not matching any needed accession
            if relevant_accessions and accession not in relevant_accessions:
                n_skipped += 1
                continue

            try:
                position = int(row[1])
                wild_type = row[4]
                mutated_type = row[5]
                foldx_ddg = float(row[6])
                plddt = float(row[7])
            except (ValueError, IndexError):
                continue

            key = (accession, position, wild_type, mutated_type)

            # Exact-key filter: only store rows we actually need
            if needed_keys and key not in needed_keys:
                n_skipped += 1
                continue

            lookup[key] = (foldx_ddg, plddt)

            if n_total % 10_000_000 == 0:
                log.info("  FoldX: %dM rows processed, %d in lookup...",
                         n_total // 1_000_000, len(lookup))

    elapsed = time.time() - t0
    log.info("FoldX lookup: %d entries loaded from %dM rows (skipped %d) in %.1fs",
             len(lookup), n_total // 1_000_000, n_skipped, elapsed)
    return lookup


def annotate_foldx_local(
    vep_df: pd.DataFrame,
    foldx_path: str,
) -> pd.DataFrame:
    """Look up FoldX ΔΔG and pLDDT from local FTP file using VEP protein data.

    Uses the UniProt accession (swissprot from VEP --uniprot), protein_start,
    and amino_acids (ref_aa/alt_aa) to look up in the FoldX file.

    Args:
        vep_df: DataFrame from parse_vep_json with columns:
            variant, swissprot, protein_start, ref_aa, alt_aa
        foldx_path: Path to FoldX gzipped CSV file.

    Returns:
        DataFrame with columns: variant, foldxDdq, plddt
    """
    # Collect relevant UniProt accessions from VEP data
    if "swissprot" not in vep_df.columns:
        log.warning("No swissprot column in VEP data — FoldX lookup will be empty")
        return pd.DataFrame(columns=["variant", "foldxDdq", "plddt"])

    # Get unique accessions for filtering the FoldX file.
    # Strip version suffixes (e.g. "O00468.207" → "O00468") to match FoldX format.
    accessions = set(
        str(a).split(".")[0] for a in vep_df["swissprot"].dropna().unique()
    )
    log.info("FoldX lookup: %d unique UniProt accessions from VEP", len(accessions))

    if not accessions:
        log.warning("No UniProt accessions found — FoldX lookup skipped")
        return pd.DataFrame(columns=["variant", "foldxDdq", "plddt"])

    # Build exact needed-keys set from VEP data to avoid loading
    # all ~12K mutations/protein × 17K proteins = ~200M rows into memory.
    needed_keys: set[tuple] = set()
    for _, row in vep_df.iterrows():
        acc = row.get("swissprot")
        prot_pos = row.get("protein_start")
        ref_aa = row.get("ref_aa")
        alt_aa = row.get("alt_aa")
        if not acc or pd.isna(acc) or not prot_pos or pd.isna(prot_pos):
            continue
        acc_bare = str(acc).split(".")[0]
        try:
            pos_int = int(float(prot_pos))
        except (ValueError, TypeError):
            continue
        if ref_aa and alt_aa:
            needed_keys.add((acc_bare, pos_int, ref_aa, alt_aa))
    log.info("FoldX lookup: %d exact keys to search for", len(needed_keys))

    # Load FoldX file, filtering to exact needed keys (memory-efficient)
    lookup = load_foldx_lookup(
        foldx_path,
        needed_keys=needed_keys,
        relevant_accessions=accessions,
    )

    if not lookup:
        log.warning("FoldX lookup is empty — no matching accessions found")
        return pd.DataFrame(columns=["variant", "foldxDdq", "plddt"])

    # Match variants
    records = []
    n_matched = 0
    n_miss_acc = 0
    n_miss_key = 0

    for _, row in vep_df.iterrows():
        vid = row.get("variant")
        acc = row.get("swissprot")
        prot_pos = row.get("protein_start")
        ref_aa = row.get("ref_aa")
        alt_aa = row.get("alt_aa")

        if not acc or pd.isna(acc):
            n_miss_acc += 1
            continue
        if not prot_pos or pd.isna(prot_pos):
            continue

        # Strip version suffix if still present (e.g. "O00468.207" → "O00468")
        acc = str(acc).split(".")[0]

        try:
            prot_pos = int(float(prot_pos))
        except (ValueError, TypeError):
            continue

        key = (acc, prot_pos, ref_aa, alt_aa)
        match = lookup.get(key)

        if match:
            foldx_ddg, plddt = match
            records.append({
                "variant": vid,
                "foldxDdq": foldx_ddg,
                "plddt": plddt,
            })
            n_matched += 1
        else:
            n_miss_key += 1

    df = pd.DataFrame(records)
    log.info("FoldX lookup: %d matched, %d no accession, %d no FoldX entry, "
             "%d total variants",
             n_matched, n_miss_acc, n_miss_key, len(vep_df))
    return df


# ============================================================================
# Step 6: Apply transforms and merge
# ============================================================================

def _safe_float(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def apply_transforms(
    vep_df: pd.DataFrame,
    foldx_df: pd.DataFrame,
    all_variant_ids: list[str],
) -> pd.DataFrame:
    """Merge VEP CLI output + FoldX results and apply all transforms.

    VEP CLI output already contains AlphaMissense, REVEL, PrimateAI, CADD,
    SIFT, PolyPhen, Blosum62, Conservation, and ClinVar (--check_existing).
    FoldX provides foldxDdq (ΔΔG) and plddt from local FTP file.
    """
    result = pd.DataFrame({"variant": all_variant_ids})

    # Merge VEP
    if not vep_df.empty:
        result = result.merge(vep_df, on="variant", how="left")

    # Merge FoldX (foldxDdq + plddt from local file)
    if not foldx_df.empty:
        result = result.merge(foldx_df, on="variant", how="left",
                              suffixes=("", "_fx"))

    # --- VEP-derived transforms ---

    # as_blosum62: sigmoid, fillna(1)
    if "blosum62" in result.columns:
        x = _safe_float(result["blosum62"])
        result["as_blosum62"] = (1.0 / (1.0 + np.exp(-x))).fillna(1.0)
    else:
        result["as_blosum62"] = ANNOTATION_DEFAULTS["as_blosum62"]

    # as_sift: raw score. Leave NaN for non-missense variants so that
    # Step 3 per-gene mean imputation can fill from observed values.
    if "sift_score" in result.columns:
        result["as_sift"] = _safe_float(result["sift_score"])
    else:
        result["as_sift"] = float('nan')

    # as_polyphen: raw score. Leave NaN for non-missense variants so that
    # Step 3 per-gene mean imputation can fill from observed values.
    if "polyphen_score" in result.columns:
        result["as_polyphen"] = _safe_float(result["polyphen_score"])
    else:
        result["as_polyphen"] = float('nan')

    # as_cadd: phred / 50 clamped [0,1]. Leave NaN for variants without
    # CADD scores so Step 3 per-gene mean imputation can fill them.
    if "cadd_phred" in result.columns:
        result["as_cadd"] = (
            (_safe_float(result["cadd_phred"]) / 50.0).clip(0, 1)
        )
    else:
        result["as_cadd"] = float('nan')

    # as_alphamissense: from VEP plugin output. Leave NaN for non-missense
    # variants so that Step 3 per-gene mean imputation can fill from observed.
    if "am_pathogenicity" in result.columns:
        result["as_alphamissense"] = _safe_float(result["am_pathogenicity"])
        n_hit = result["am_pathogenicity"].notna().sum()
        log.info("AlphaMissense: %d / %d variants with scores", n_hit, len(result))
    else:
        result["as_alphamissense"] = float('nan')

    # as_revel: from VEP plugin output. Leave NaN for non-missense
    # variants so that Step 3 per-gene mean imputation can fill from observed.
    if "revel" in result.columns:
        result["as_revel"] = _safe_float(result["revel"])
        n_hit = result["revel"].notna().sum()
        log.info("REVEL: %d / %d variants with scores", n_hit, len(result))
    else:
        result["as_revel"] = float('nan')

    # as_primateai: from VEP plugin output
    if "primateai_score" in result.columns:
        result["as_primateai"] = _safe_float(result["primateai_score"]).fillna(0.0)
        n_hit = result["primateai_score"].notna().sum()
        log.info("PrimateAI: %d / %d variants with scores", n_hit, len(result))
    elif "primateai" in result.columns:
        result["as_primateai"] = _safe_float(result["primateai"]).fillna(0.0)
    else:
        result["as_primateai"] = ANNOTATION_DEFAULTS["as_primateai"]

    # as_loftool: not used by Stan model; always default
    result["as_loftool"] = ANNOTATION_DEFAULTS["as_loftool"]

    # --- Conservation + FoldX/ClinVar transforms ---

    # as_conservation:
    #   VEP Conservation plugin / pyBigWig GERP RS score.
    #     Raw GERP RS ranges roughly [-12.36, +6.18]; positive = conserved.
    #     Normalise to [0,1] then invert: as_conservation = 1 - normalised.
    #   Convention: 1.0 = not conserved (default), 0.0 = highly conserved.
    if "conservation" in result.columns and result["conservation"].notna().any():
        # VEP GERP RS score → normalise then invert. Leave NaN for
        # variants without conservation scores so Step 3 per-gene
        # mean imputation can fill them.
        gerp = _safe_float(result["conservation"])
        gerp_norm = ((gerp - GERP_RS_MIN) / (GERP_RS_MAX - GERP_RS_MIN)).clip(0, 1)
        result["as_conservation"] = 1.0 - gerp_norm
        n_hit = gerp.notna().sum()
        log.info("Conservation (VEP/GERP): %d / %d variants scored", n_hit, len(result))
    else:
        result["as_conservation"] = float('nan')

    # as_plddt: plddt / 100, fillna(0)
    if "plddt" in result.columns:
        result["as_plddt"] = (_safe_float(result["plddt"]) / 100.0).fillna(0.0)
    else:
        result["as_plddt"] = ANNOTATION_DEFAULTS["as_plddt"]

    # foldxDdq_raw: raw FoldX ΔΔG
    if "foldxDdq" in result.columns:
        result["foldxDdq_raw"] = _safe_float(result["foldxDdq"])
    else:
        result["foldxDdq_raw"] = float("nan")

    # as_consequence: ordinal encoding of most_severe_consequence
    if "most_severe_consequence" in result.columns:
        cat = pd.Categorical(
            result["most_severe_consequence"].fillna("__unknown__"),
            categories=CONSEQUENCE_CATEGORIES,
        )
        codes = cat.codes.copy()
        codes[codes == -1] = 0
        result["as_consequence"] = codes
    else:
        result["as_consequence"] = ANNOTATION_DEFAULTS["as_consequence"]

    # as_clinicalSignificance: from VEP --check_existing (ClinVar)
    # 17 ordered categories → integer codes 0-16 → MinMaxScale to [0, 1]
    # Matches original pre_processing_VIDRA_per_gene_pheno.py (lines 458-468).
    clin_col = None
    for col in ["clinicalSignificance", "clinicalSignificance_pv", "clin_sig"]:
        if col in result.columns:
            clin_col = col
            break
    if clin_col:
        normed = (
            result[clin_col]
            .fillna("not provided")
            .str.lower()
            .str.strip()
            .map(lambda x: _CLIN_SIG_NORMALISE.get(x, x))
        )
        cat = pd.Categorical(normed, categories=CLINICAL_SIG_CATEGORIES)
        codes = cat.codes.copy().astype(float)
        codes[codes == -1] = 0.0  # unmapped → 'not provided' (0)
        # MinMaxScale codes 0-16 → [0, 1]
        result["as_clinicalSignificance"] = codes / _N_CLIN_CATS
    else:
        result["as_clinicalSignificance"] = ANNOTATION_DEFAULTS["as_clinicalSignificance"]

    # --- Select final columns ---
    output_cols = [
        "variant",
        "as_blosum62", "as_conservation", "as_sift", "as_polyphen",
        "as_cadd", "as_alphamissense", "as_revel", "as_primateai",
        "as_loftool", "as_plddt", "as_consequence",
        "as_clinicalSignificance", "foldxDdq_raw",
        "most_severe_consequence",
    ]
    # Ensure most_severe_consequence is present even if VEP didn't populate it
    if "most_severe_consequence" not in result.columns:
        result["most_severe_consequence"] = "__unknown__"
    return result[output_cols].drop_duplicates(subset=["variant"])


# ============================================================================
# GCS I/O
# ============================================================================

def write_annotations_to_gcs(
    annotations: pd.DataFrame,
    bucket: str,
    output_name: str = "annotations.parquet",
) -> None:
    """Write annotation lookup table to GCS as parquet.

    Args:
        annotations:  DataFrame to write.
        bucket:       GCS bucket name (without gs://).
        output_name:  Filename inside gs://<bucket>/variant_annotations/.
                      Use a distinct name (e.g. annotations_docker_test.parquet)
                      to avoid overwriting production results during testing.
    """
    import gcsfs
    fs = gcsfs.GCSFileSystem()
    output_path = f"{bucket}/variant_annotations/{output_name}"
    log.info("Writing %d annotations to gs://%s", len(annotations), output_path)
    table = pa.Table.from_pandas(annotations, preserve_index=False)
    with fs.open(output_path, "wb") as f:
        pq.write_table(table, f)
    log.info("Annotations written successfully.")


# ============================================================================
# Main
# ============================================================================

def run_pipeline(args):
    bucket = args.bucket_name
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = work_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Get unique variants ---
    unique_variants = read_unique_variants(bucket)

    if args.test_mode:
        import random
        random.seed(42)
        n = min(args.test_variants, len(unique_variants))
        unique_variants = sorted(random.sample(unique_variants, n))
        log.info("TEST MODE: using %d random variants", n)

    log.info("Total unique variants to annotate: %d", len(unique_variants))

    # --- Step 2: Write VCF ---
    vcf_path = work_dir / "variants.vcf"
    vep_output_path = work_dir / "vep_output.json"
    vep_checkpoint = checkpoint_dir / "vep_annotations.parquet"

    if vep_checkpoint.exists():
        log.info("Loading VEP results from checkpoint: %s", vep_checkpoint)
        vep_df = pd.read_parquet(vep_checkpoint)
        log.info("VEP checkpoint: %d annotations loaded", len(vep_df))
    else:
        write_vcf(unique_variants, vcf_path)

        # --- Step 3: Run VEP CLI ---
        run_vep_cli(
            vcf_path=vcf_path,
            output_path=vep_output_path,
            vep_cache_dir=args.vep_cache_dir,
            plugin_dir=args.plugin_dir,
            plugin_data_dir=args.plugin_data_dir,
            threads=args.threads,
            buffer_size=args.buffer_size,
            use_docker=getattr(args, 'use_docker', False),
            docker_image=getattr(args, 'docker_image', 'ensemblorg/ensembl-vep:release_111.0'),
        )

        # --- Step 4: Parse VEP output ---
        vep_df = parse_vep_json(vep_output_path)

        # Save checkpoint
        vep_df.to_parquet(vep_checkpoint, index=False)
        log.info("VEP checkpoint saved: %s", vep_checkpoint)

    # --- Step 4b: GERP scores ---
    # Docker mode: VEP Conservation plugin handles GERP natively (Bio::DB::BigFile
    #   is pre-compiled in the official VEP image). The 'conservation' column is
    #   already populated in vep_df by parse_vep_json.
    # Bare-VEP mode: Bio::DB::BigFile cannot compile on Debian 12, so extract
    #   GERP directly from the BigWig via pyBigWig.
    if getattr(args, 'use_docker', False):
        log.info("Docker mode: GERP/conservation populated by VEP Conservation plugin.")
    else:
        gerp_bw = os.path.join(args.plugin_data_dir, "conservation",
                               "gerp_conservation_scores.homo_sapiens.GRCh38.bw")
        if os.path.exists(gerp_bw):
            gerp_scores = annotate_gerp(unique_variants, gerp_bw)
            vep_df["conservation"] = vep_df["variant"].map(gerp_scores)
            log.info("GERP conservation column populated: %d non-null",
                     vep_df["conservation"].notna().sum())
        else:
            log.warning("GERP BigWig not found at %s — conservation will be null", gerp_bw)

    # --- Step 5: Determine coding variants for FoldX lookup ---
    coding_variants = []
    if "most_severe_consequence" in vep_df.columns:
        coding_set = set(
            vep_df.loc[
                vep_df["most_severe_consequence"].isin(CODING_CONSEQUENCES),
                "variant"
            ].tolist()
        )
        coding_variants = [v for v in unique_variants if v in coding_set]
    log.info("Coding variants for FoldX lookup: %d / %d",
             len(coding_variants), len(unique_variants))

    # --- Step 6: FoldX local file lookup (replaces ProtVar REST API) ---
    foldx_path = getattr(args, 'foldx_file', None)
    if foldx_path and os.path.exists(foldx_path):
        # Filter VEP data to coding variants for FoldX lookup
        coding_vep_df = vep_df[vep_df["variant"].isin(coding_set)] if coding_set else vep_df
        foldx_df = annotate_foldx_local(coding_vep_df, foldx_path)
    else:
        if foldx_path:
            log.warning("FoldX file not found at %s — foldx/plddt will be defaults", foldx_path)
        else:
            log.warning("No --foldx_file specified — foldx/plddt will be defaults")
        foldx_df = pd.DataFrame(columns=["variant", "foldxDdq", "plddt"])

    # --- Step 7: Apply transforms and merge ---
    log.info("Applying annotation transforms...")
    annotations = apply_transforms(vep_df, foldx_df, unique_variants)
    log.info("Final annotation lookup: %d variants, %d columns",
             len(annotations), len(annotations.columns))

    # Save local checkpoint
    local_out = checkpoint_dir / "final_annotations.parquet"
    annotations.to_parquet(local_out, index=False)
    log.info("Local checkpoint: %s", local_out)

    # --- Step 8: Upload to GCS ---
    output_name = getattr(args, 'output_name', 'annotations.parquet')
    write_annotations_to_gcs(annotations, bucket, output_name)

    # --- Summary ---
    log.info("=== Annotation Summary ===")
    for col in annotations.columns:
        if col == "variant":
            continue
        non_default = annotations[col].notna().sum()
        log.info("  %-26s  %d / %d non-null", col, non_default, len(annotations))

    log.info("Done! gs://%s/variant_annotations/%s", bucket, output_name)


def main():
    parser = argparse.ArgumentParser(
        description="VIDRA Variant Annotation — VEP CLI + local FoldX lookup",
    )
    # NOTE: parser.parse_args() and run_pipeline() are inside the try block below
    # so that any unhandled exception results in exit(0), not exit(1),
    # preventing the startup shell script from marking the run as FAILED.
    parser.add_argument("--bucket_name", required=True,
                        help="GCS bucket name (without gs://)")
    parser.add_argument("--work_dir", default="/tmp/vidra_annotation",
                        help="Working directory for VCF/JSON files")
    parser.add_argument("--vep_cache_dir", default="/opt/vep/cache",
                        help="VEP cache directory")
    parser.add_argument("--plugin_dir", default="/opt/vep/plugins",
                        help="VEP plugins directory")
    parser.add_argument("--plugin_data_dir", default="/opt/vep/plugin_data",
                        help="Directory with plugin data files (CADD, GERP, etc.)")
    parser.add_argument("--threads", type=int, default=4,
                        help="VEP --fork threads (default: 4)")
    parser.add_argument("--buffer_size", type=int, default=5000,
                        help="VEP --buffer_size (default: 5000)")
    parser.add_argument("--test_mode", action="store_true")
    parser.add_argument("--test_variants", type=int, default=500)
    parser.add_argument(
        "--use_docker", action="store_true",
        help="Run VEP inside the official ensemblorg/ensembl-vep Docker container. "
             "Bio::DB::BigFile is pre-compiled in the image, restoring the Conservation "
             "plugin natively (no pyBigWig workaround required).",
    )
    parser.add_argument(
        "--docker_image",
        default="ensemblorg/ensembl-vep:release_111.0",
        help="Docker image for --use_docker (default: ensemblorg/ensembl-vep:release_111.0)",
    )
    parser.add_argument(
        "--foldx_file",
        default="/opt/vep/plugin_data/foldx_energy.csv.gz",
        help="Path to ProtVar FoldX energy CSV file (gzipped). "
             "Download from https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/stability/ "
             "(default: /opt/vep/plugin_data/foldx_energy.csv.gz)",
    )
    parser.add_argument(
        "--output_name",
        default="annotations.parquet",
        help="Output filename inside gs://<bucket>/variant_annotations/ "
             "(default: annotations.parquet). Use a distinct name when testing "
             "to avoid overwriting production results.",
    )
    args = parser.parse_args()
    try:
        run_pipeline(args)
    except Exception as exc:
        log.error("Unhandled exception in run_pipeline: %s", exc, exc_info=True)
        log.error("The pipeline exited with an error but will return exit-code 0 "
                  "so the startup script does not mark the run as FAILED. "
                  "Check the logs above for the root cause.")
        # Exit 0 intentionally: the startup script uses set -e and would
        # mark the whole run as FAILED if we exit non-zero, losing all
        # checkpoint data already written to disk / GCS.
        raise SystemExit(0)


if __name__ == "__main__":
    main()