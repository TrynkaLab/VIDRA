#!/usr/bin/env bash
# =============================================================================
# run_annotation_vm.sh — Create a GCE VM that autonomously runs VEP + FoldX
# =============================================================================
#
# Fully autonomous: the VM runs a startup script that installs VEP v111,
# downloads all plugin data + FoldX energy file, executes
# annotate_variants_cli.py, uploads results to GCS, and optionally
# self-deletes. No SSH session required.
#
# Progress monitoring:
#   gsutil cat gs://vidra-2-0/annotation_status/status.txt
#   gcloud compute instances get-serial-port-output vidra-annotation-vm \
#     --project=open-targets-genetics-dev --zone=europe-west1-b | tail -50
#
# Usage:
#   bash scripts/pyspark_scripts/run_annotation_vm.sh              # full run + auto-delete
#   bash scripts/pyspark_scripts/run_annotation_vm.sh --no-delete  # keep VM for debugging
#
# Estimated cost: ~$2-5 (e2-standard-8, SSD, ~1.5-3 hours total)
# Estimated time breakdown:
#   - System setup + VEP install:  ~10 min
#   - Plugin data downloads:       ~30-90 min (CADD ~80GB is the bottleneck)
#   - VEP CLI annotation:          ~30-60 min (1.6M variants, 4 parallel chunks × 2 forks)
#   - FoldX local lookup:          ~5 min (1.2M coding variants)
# =============================================================================

set -euo pipefail

# --- Configuration ---
PROJECT="open-targets-genetics-dev"
ZONE="europe-west1-b"
VM_NAME="vidra-annotation-vm"
MACHINE_TYPE="e2-standard-8"     # 8 vCPU, 32 GB RAM
BOOT_DISK_SIZE="250GB"           # VEP cache ~15GB, CADD ~85GB, working space
BUCKET="vidra-2-0"
VEP_VERSION="111"
AUTO_DELETE="true"
OUTPUT_SUFFIX=""

for arg in "$@"; do
    case "$arg" in
        --no-delete)    AUTO_DELETE="false" ;;
        --suffix=*)     OUTPUT_SUFFIX="${arg#*=}" ;;
    esac
done

echo "=== VIDRA Annotation VM ==="
echo "  Project:      $PROJECT"
echo "  Zone:         $ZONE"
echo "  VM:           $VM_NAME"
echo "  Machine:      $MACHINE_TYPE"
echo "  Disk:         $BOOT_DISK_SIZE"
echo "  VEP version:  $VEP_VERSION"
echo "  Auto-delete:  $AUTO_DELETE"
echo "  Suffix:       ${OUTPUT_SUFFIX:-(none)}"
echo ""

# --- Step 0: Upload scripts to GCS ---
echo "Uploading scripts to GCS..."
gsutil cp scripts/pyspark_scripts/annotate_variants_cli.py \
    "gs://${BUCKET}/scripts/annotate_variants_cli.py"

# --- Step 1: Create the startup script as a temp file ---
# Using --metadata-from-file to avoid shell quoting issues with --metadata=
STARTUP_FILE=$(mktemp /tmp/vidra_startup_XXXXXX)
cat > "$STARTUP_FILE" <<'STARTUP_EOF'
#!/bin/bash
set -euo pipefail

BUCKET="vidra-2-0"
OUTPUT_SUFFIX_FLAG="__OUTPUT_SUFFIX_PLACEHOLDER__"
LOG="/var/log/vidra-annotation.log"
STATUS_GCS="gs://${BUCKET}/annotation_status/status.txt"
PLUGIN_DATA="/opt/vep/plugin_data"
VEP_HOME="/opt/vep"

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"
}

update_status() {
    local msg="$1"
    log_msg "STATUS: $msg"
    echo "$msg — $(date -u '+%Y-%m-%d %H:%M:%S UTC')" | gcloud storage cp - "$STATUS_GCS"
}

on_error() {
    update_status "FAILED at line $1 — check serial console or $LOG"
    exit 1
}
trap 'on_error $LINENO' ERR

# -------------------------------------------------------------------
# Phase 1: System packages
# -------------------------------------------------------------------
update_status "Phase 1/6: Installing system packages"
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq \
    docker.io \
    tabix samtools unzip wget curl \
    python3 python3-pip python3-venv git

log_msg "System packages installed"

# -------------------------------------------------------------------
# Phase 2: Pull VEP v111 Docker image + download GRCh38 cache
# -------------------------------------------------------------------
update_status "Phase 2/6: Pulling VEP v111 Docker image + downloading GRCh38 cache"
mkdir -p "${VEP_HOME}/cache"

# Enable and start Docker daemon (installed in Phase 1)
systemctl start docker
systemctl enable docker
log_msg "Docker version: $(docker --version)"

# Pull the official VEP v111 image (~5 GB download).
log_msg "Pulling ensemblorg/ensembl-vep:release_111.0..."
docker pull ensemblorg/ensembl-vep:release_111.0 2>&1 | tee -a "$LOG"
log_msg "VEP Docker image pulled"

# Quick sanity check — confirm vep binary responds inside the container
docker run --rm ensemblorg/ensembl-vep:release_111.0 vep --help 2>&1 | head -3 | tee -a "$LOG"

# Download VEP GRCh38 cache directly from Ensembl FTP (~15 GB, ~15-20 min).
log_msg "Downloading VEP GRCh38 cache (this takes ~15-20 min)..."
mkdir -p "${VEP_HOME}/cache"
wget -q --show-progress \
    -O "${VEP_HOME}/cache/homo_sapiens_vep_111_GRCh38.tar.gz" \
    "https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz" \
    2>&1 | tee -a "$LOG"
log_msg "Extracting VEP cache (~5 min)..."
tar -xzf "${VEP_HOME}/cache/homo_sapiens_vep_111_GRCh38.tar.gz" -C "${VEP_HOME}/cache/" 2>&1 | tee -a "$LOG"
rm -f "${VEP_HOME}/cache/homo_sapiens_vep_111_GRCh38.tar.gz"
log_msg "VEP GRCh38 cache ready at ${VEP_HOME}/cache"
ls "${VEP_HOME}/cache/" | tee -a "$LOG"

# -------------------------------------------------------------------
# Phase 3: Copy plugin data files from GCS
# -------------------------------------------------------------------
# All plugin data is pre-staged on gs://vidra-2-0/raw_data/
# This is MUCH faster than downloading from external servers.
update_status "Phase 3/6: Copying plugin data from GCS to local disk"
mkdir -p "$PLUGIN_DATA"/{CADD,conservation}

# -- AlphaMissense (599 MB) ---
log_msg "Copying AlphaMissense..."
gcloud storage cp "gs://${BUCKET}/raw_data/AlphaMissense_hg38.tsv.gz" "$PLUGIN_DATA/"
tabix -s 1 -b 2 -e 2 -f "$PLUGIN_DATA/AlphaMissense_hg38.tsv.gz" || true
log_msg "AlphaMissense ready"

# -- REVEL GRCh38 (filtered grch38_pos only, sorted + indexed on col3) ---
log_msg "Copying REVEL (GRCh38 re-indexed)..."
gcloud storage cp "gs://${BUCKET}/raw_data/new_tabbed_revel_grch38.tsv.gz" "$PLUGIN_DATA/"
gcloud storage cp "gs://${BUCKET}/raw_data/new_tabbed_revel_grch38.tsv.gz.tbi" "$PLUGIN_DATA/"
log_msg "REVEL ready"

# -- PrimateAI (839 MB) ---
log_msg "Copying PrimateAI..."
gcloud storage cp "gs://${BUCKET}/raw_data/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz" "$PLUGIN_DATA/"
tabix -s 1 -b 2 -e 2 -f "$PLUGIN_DATA/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz" || true
log_msg "PrimateAI ready"

# -- CADD (~80 GB) ---
# gcloud storage uses gRPC and parallel composite downloads by default —
# much faster than gsutil cp for large single-object transfers.
log_msg "Copying CADD whole_genome_SNVs.tsv.gz (~80GB from GCS)..."
gcloud storage cp "gs://${BUCKET}/raw_data/CADD/whole_genome_SNVs.tsv.gz" "$PLUGIN_DATA/CADD/"
gcloud storage cp "gs://${BUCKET}/raw_data/CADD/whole_genome_SNVs.tsv.gz.tbi" "$PLUGIN_DATA/CADD/"
log_msg "CADD SNVs copied"

# -- Conservation / GERP (8.9 GB) ---
log_msg "Copying GERP conservation scores..."
gcloud storage cp "gs://${BUCKET}/raw_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
    "$PLUGIN_DATA/conservation/"
log_msg "GERP conservation ready"

# -- FoldX energy file from ProtVar FTP (1.5 GB) ---
# Replaces the unreliable ProtVar REST API for foldxDdq + plddt.
# Source: https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/stability/
log_msg "Downloading FoldX energy file from ProtVar FTP (~1.5GB)..."
if [[ -f "$PLUGIN_DATA/foldx_energy.csv.gz" ]]; then
    log_msg "FoldX file already exists, skipping download"
else
    # Try GCS copy first (if pre-staged), fall back to FTP download
    if gcloud storage cp "gs://${BUCKET}/raw_data/foldx_energy.csv.gz" "$PLUGIN_DATA/" 2>/dev/null; then
        log_msg "FoldX copied from GCS"
    else
        log_msg "FoldX not on GCS — downloading from EBI FTP..."
        wget -q --show-progress \
            -O "$PLUGIN_DATA/foldx_energy.csv.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/stability/2025.02.10_foldx_energy.csv.gz" \
            2>&1 | tee -a "$LOG"
        log_msg "FoldX downloaded from FTP"
        # Stage on GCS for future runs
        gcloud storage cp "$PLUGIN_DATA/foldx_energy.csv.gz" "gs://${BUCKET}/raw_data/foldx_energy.csv.gz" || true
    fi
fi
log_msg "FoldX energy file ready"

log_msg "Plugin data summary:"
ls -lhR "$PLUGIN_DATA/" 2>&1 | tee -a "$LOG"

# -------------------------------------------------------------------
# Phase 4: Python environment
# -------------------------------------------------------------------
update_status "Phase 4/6: Setting up Python environment"
python3 -m venv "$VEP_HOME/pyenv"
source "$VEP_HOME/pyenv/bin/activate"
pip install --quiet --upgrade pip
pip install --quiet pandas numpy pyarrow gcsfs
log_msg "Python venv ready"

# -------------------------------------------------------------------
# Phase 5: Download annotation script
# -------------------------------------------------------------------
update_status "Phase 5/6: Downloading annotation script from GCS"
gsutil cp "gs://${BUCKET}/scripts/annotate_variants_cli.py" \
    "$VEP_HOME/annotate_variants_cli.py"
log_msg "Script downloaded"

# -------------------------------------------------------------------
# Phase 6: Run VEP CLI + FoldX local annotation
# -------------------------------------------------------------------
update_status "Phase 6/6: Running annotation (VEP CLI + FoldX local lookup)"
source "$VEP_HOME/pyenv/bin/activate"

python3 "$VEP_HOME/annotate_variants_cli.py" \
    --bucket_name "$BUCKET" \
    --work_dir /tmp/vidra_annotation \
    --vep_cache_dir "$VEP_HOME/cache" \
    --plugin_dir "$VEP_HOME/plugins" \
    --plugin_data_dir "$PLUGIN_DATA" \
    --threads 8 \
    --buffer_size 10000 \
    --use_docker \
    --docker_image ensemblorg/ensembl-vep:release_111.0 \
    --foldx_file "$PLUGIN_DATA/foldx_energy.csv.gz" \
    --vep_parallel 4 \
    --output_name variant_annotations.parquet \
    --output_suffix "$OUTPUT_SUFFIX_FLAG" \
    2>&1 | tee -a "$LOG"

log_msg "Annotation script completed"

# -------------------------------------------------------------------
# Upload log and clean up
# -------------------------------------------------------------------
update_status "Uploading logs and finishing"
gsutil cp "$LOG" "gs://${BUCKET}/annotation_status/annotation_run.log"
log_msg "Log uploaded to gs://${BUCKET}/annotation_status/annotation_run.log"

update_status "COMPLETED — annotations at gs://${BUCKET}/variant_annotations${OUTPUT_SUFFIX_FLAG}/"

# --- Self-delete if configured ---
AUTO_DELETE_FLAG="__AUTO_DELETE_PLACEHOLDER__"
if [[ "$AUTO_DELETE_FLAG" == "true" ]]; then
    log_msg "Auto-deleting VM in 60 seconds..."
    sleep 60
    VM_NAME=$(curl -s -H "Metadata-Flavor: Google" \
        http://metadata.google.internal/computeMetadata/v1/instance/name)
    VM_ZONE=$(curl -s -H "Metadata-Flavor: Google" \
        http://metadata.google.internal/computeMetadata/v1/instance/zone | awk -F/ '{print $NF}')
    gcloud compute instances delete "$VM_NAME" --zone="$VM_ZONE" --quiet || true
fi
STARTUP_EOF

# Replace the auto-delete placeholder with the actual setting
sed -i '' "s/__AUTO_DELETE_PLACEHOLDER__/$AUTO_DELETE/g" "$STARTUP_FILE"
sed -i '' "s/__OUTPUT_SUFFIX_PLACEHOLDER__/$OUTPUT_SUFFIX/g" "$STARTUP_FILE"

# --- Step 2: Create the VM with the startup script ---
echo "Creating VM: $VM_NAME ..."
gcloud compute instances create "$VM_NAME" \
    --project="$PROJECT" \
    --zone="$ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --boot-disk-size="$BOOT_DISK_SIZE" \
    --boot-disk-type=pd-ssd \
    --image-family=debian-12 \
    --image-project=debian-cloud \
    --scopes=cloud-platform \
    --metadata-from-file=startup-script="$STARTUP_FILE" \
    --quiet

rm -f "$STARTUP_FILE"

echo ""
echo "=== VM created and running autonomously ==="
echo ""
echo "The VM is now executing the full annotation pipeline."
echo "No SSH connection needed — it runs entirely via startup script."
echo ""
echo "Monitor progress:"
echo "  gsutil cat gs://${BUCKET}/annotation_status/status.txt"
echo ""
echo "View full logs:"
echo "  gsutil cat gs://${BUCKET}/annotation_status/annotation_run.log"
echo ""
echo "Serial console (live):"
echo "  gcloud compute instances get-serial-port-output $VM_NAME \\"
echo "    --project=$PROJECT --zone=$ZONE | tail -100"
echo ""
echo "SSH in if needed:"
echo "  gcloud compute ssh $VM_NAME --project=$PROJECT --zone=$ZONE"
echo ""
if [[ "$AUTO_DELETE" == "true" ]]; then
    echo "VM will self-delete after completion."
else
    echo "VM will remain running. Delete manually when done:"
    echo "  gcloud compute instances delete $VM_NAME --project=$PROJECT --zone=$ZONE --quiet"
fi
