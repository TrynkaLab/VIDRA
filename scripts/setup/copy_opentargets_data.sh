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

# 5. Burden (Evidence) - gene burden tests from OT Platform
# Note: Source updated to 22.04 per user request
echo "Copying Burden Data (approx 10MB)..."
gcloud storage rsync \
  "gs://open-targets-data-releases/22.04/output/etl/parquet/evidence/sourceId=gene_burden" \
  "$DEST/burden" \
  --recursive \
  --delete-unmatched-destination-objects \
  --billing-project=$PROJECT

# 6. ClinVar (Evidence) - ClinVar variant evidence from OT Platform
# Note: Source updated to 24.03 per user request
echo "Copying ClinVar Data (approx 200MB)..."
gcloud storage rsync \
  "gs://open-targets-data-releases/24.03/output/etl/parquet/evidence/sourceId=eva" \
  "$DEST/clinvar" \
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
echo "  $DEST/burden/   - gene burden evidence"
echo "  $DEST/clinvar/  - ClinVar evidence (eva)"
