#!/bin/bash
set -e

# Copy Open Targets data releases into the VIDRA bucket for PySpark processing.
# Source: OT Genetics 22.09 + OT Platform 22.09
# Target: gs://vidra-2-0/raw_data/open_targets/
#
# Prerequisites:
#   - gcloud auth login
#   - Access to open-targets-genetics-dev project (for requester-pays buckets)
#
# Usage:
#   bash scripts/copy_opentargets_data.sh

# Target Bucket
BUCKET="gs://vidra-2-0"
DEST="$BUCKET/raw_data/open_targets"
PROJECT="open-targets-genetics-dev"

echo "Starting data copy to $DEST..."

# 1. Coloc (Genetics) - variant-disease colocalisation
echo "Copying Coloc Data (approx 600MB)..."
gcloud storage rsync \
  "gs://open-targets-genetics-data-releases/open-targets-genetics-releases/22.09/v2d_coloc" \
  "$DEST/coloc" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 2. GWAS (Genetics) - summary association statistics
echo "Copying GWAS Data (approx 25GB)..."
gcloud storage rsync \
  "gs://open-targets-genetics-data-releases/open-targets-genetics-releases/22.09/sa/gwas" \
  "$DEST/gwas" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 3. Variants (Genetics) - variant index with gene assignments & consequences
echo "Copying Variants Data (approx 8GB)..."
gcloud storage rsync \
  "gs://open-targets-genetics-data-releases/open-targets-genetics-releases/22.09/variant-index" \
  "$DEST/variants" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 4. Studies (Genetics) - study index with trait_efos for phenotype mapping
echo "Copying Studies Data (approx 5MB)..."
gcloud storage rsync \
  "gs://open-targets-genetics-data-releases/open-targets-genetics-releases/22.09/lut/study-index" \
  "$DEST/studies" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 5. AZ PheWAS Burden (raw collapsing model data, replaces OT gene_burden)
# Binary and quantitative parquets have different schemas (BinOddsRatio vs beta).
# These must be uploaded from local downloads — not available on GCS.
# Upload with:
#   gcloud storage cp --recursive /path/to/azphewas-com-450k-phewas-binary/ \
#     gs://vidra-2-0/raw_data/az_phewas_burden/binary/
#   gcloud storage cp --recursive /path/to/azphewas-com-450k-phewas-quantitative/ \
#     gs://vidra-2-0/raw_data/az_phewas_burden/quantitative/
echo "Checking AZ PheWAS Burden Data..."
AZ_BURDEN="$BUCKET/raw_data/az_phewas_burden"
if gcloud storage ls "$AZ_BURDEN/binary/" > /dev/null 2>&1 && \
   gcloud storage ls "$AZ_BURDEN/quantitative/" > /dev/null 2>&1; then
  echo "  AZ PheWAS burden data found at $AZ_BURDEN"
else
  echo "  WARNING: AZ PheWAS burden data not found at $AZ_BURDEN"
  echo "  Upload binary/ and quantitative/ parquet directories manually."
fi

# 6. ClinVar (Evidence) - ClinVar variant evidence from OT Platform
# Note: Source updated to 24.03 per user request
echo "Copying ClinVar Data (approx 200MB)..."
gcloud storage rsync \
  "gs://open-targets-data-releases/24.03/output/etl/parquet/evidence/sourceId=eva" \
  "$DEST/clinvar" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 7. Molecule (OT Platform) - molecule index for drug target mapping
echo "Copying Molecule Data..."
gcloud storage rsync \
  "gs://open-targets-data-releases/24.03/output/etl/parquet/molecule" \
  "$DEST/molecule" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 8. Mechanism of Action (OT Platform) - mechanism of action evidence for drug targets
echo "Copying Mechanism of Action Data..."
gcloud storage rsync \
  "gs://open-targets-data-releases/24.03/output/etl/parquet/mechanismOfAction" \
  "$DEST/mechanismOfAction" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 9. Targets (OT Platform) - target index for symbol -> ENSG mapping
echo "Copying Targets Data..."
gcloud storage rsync \
  "gs://open-targets-data-releases/24.03/output/etl/parquet/targets" \
  "$DEST/targets" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 10. Diseases (OT Platform) - disease index for obsolete EFO ID remapping
echo "Copying Diseases Data..."
gcloud storage rsync \
  "gs://open-targets-data-releases/24.03/output/etl/parquet/diseases" \
  "$DEST/diseases" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

echo "All transfers complete!"
echo ""
echo "Expected layout:"
echo "  $DEST/coloc/    - colocalisation (v2d_coloc)"
echo "  $DEST/gwas/     - GWAS summary stats (sa_gwas)"
echo "  $DEST/variants/ - variant index (gene_id_prot_coding, consequences)"
echo "  $DEST/studies/  - study index (study_id -> trait_efos)"
echo "  $BUCKET/raw_data/az_phewas_burden/ - AZ PheWAS collapsing model (binary + quantitative)"
echo "  $DEST/clinvar/  - ClinVar evidence (eva)"
echo "  $DEST/molecule/ - molecule index (drug targets)"
echo "  $DEST/mechanismOfAction/ - mechanism of action evidence"
echo "  $DEST/targets/  - target index (symbol -> ENSG mapping)"
echo "  $DEST/diseases/ - disease index (obsolete EFO ID remapping)"