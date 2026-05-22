#!/bin/bash
# ==============================================================================
# TRUST4 Batch Analysis Script for HAE Multi-omics Study
# ==============================================================================
# This script runs TRUST4 on all BAM files to extract TCR/BCR clonotypes
# Prerequisites:
#   1. TRUST4 installed in multiomics conda environment
#   2. Reference files downloaded (run download_trust4_references.sh first)
#   3. External HDD mounted with BAM files
# ==============================================================================

set -e

# Configuration
BASE_DIR="/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
BAM_DIR="/Volumes/阿热-实验数据备份/工作/1.肝包虫/1.泡型包虫/1.2024-11-泡球蚴/1.转录组学测序/02.Bam"
REF_DIR="${BASE_DIR}/analysis/data/reference/trust4"
OUTPUT_DIR="${BASE_DIR}/analysis/results/enhancement_tcr_bcr/trust4_output"
LOG_DIR="${BASE_DIR}/analysis/results/enhancement_tcr_bcr/logs"

# Reference files
REF_GENOME="${REF_DIR}/hg38_bcrtcr.fa"
REF_IMGT="${REF_DIR}/human_IMGT+C.fa"

# TRUST4 path (use run-trust4 wrapper which supports --ref and --od)
TRUST4="/Users/rishat/miniforge3/envs/multiomics/bin/run-trust4"
TRUST4_SIMPLEREP="/Users/rishat/miniforge3/envs/multiomics/bin/trust-simplerep.pl"

# Number of threads
THREADS=4

echo "=== TRUST4 Batch Analysis for HAE Multi-omics Study ==="
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Check prerequisites
echo "=== Checking Prerequisites ==="

# 1. Check TRUST4
if [ ! -f "${TRUST4}" ]; then
    echo "ERROR: TRUST4 not found at ${TRUST4}"
    echo "Please install TRUST4: conda install -n multiomics bioconda::trust4"
    exit 1
fi
echo "✓ TRUST4 found: ${TRUST4}"

# 2. Check reference files
if [ ! -f "${REF_GENOME}" ]; then
    echo "ERROR: Reference genome not found: ${REF_GENOME}"
    echo "Please run: bash ${BASE_DIR}/analysis/scripts/download_trust4_references.sh"
    exit 1
fi
echo "✓ Reference genome found: ${REF_GENOME}"

if [ ! -f "${REF_IMGT}" ]; then
    echo "ERROR: IMGT reference not found: ${REF_IMGT}"
    echo "Please run: bash ${BASE_DIR}/analysis/scripts/download_trust4_references.sh"
    exit 1
fi
echo "✓ IMGT reference found: ${REF_IMGT}"

# 3. Check BAM directory
if [ ! -d "${BAM_DIR}" ]; then
    echo "ERROR: BAM directory not found: ${BAM_DIR}"
    echo "Please mount the external HDD: 阿热-实验数据备份"
    exit 1
fi
echo "✓ BAM directory found: ${BAM_DIR}"

# Count BAM files
BAM_COUNT=$(ls -1 "${BAM_DIR}"/*.bam 2>/dev/null | wc -l | tr -d ' ')
echo "✓ Found ${BAM_COUNT} BAM files"

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
echo "✓ Output directories created"
echo ""

# Get list of BAM files
BAM_FILES=$(ls -1 "${BAM_DIR}"/*.bam 2>/dev/null)

# Process each BAM file
echo "=== Running TRUST4 Analysis ==="
echo "Processing ${BAM_COUNT} samples..."
echo ""

processed=0
skipped=0
failed=0

for bam_file in ${BAM_FILES}; do
    sample_name=$(basename "${bam_file}" .bam)
    sample_output="${OUTPUT_DIR}/${sample_name}"
    log_file="${LOG_DIR}/${sample_name}.log"
    
    # Check if already processed
    if [ -f "${sample_output}_report.tsv" ]; then
        echo "[SKIP] ${sample_name} - already processed"
        ((skipped++))
        continue
    fi
    
    echo "[RUN] ${sample_name}"
    echo "  BAM: ${bam_file}"
    echo "  Output: ${sample_output}"
    
    # Run TRUST4
    if "${TRUST4}" \
        -b "${bam_file}" \
        -f "${REF_GENOME}" \
        --ref "${REF_IMGT}" \
        -o "${sample_name}" \
        --od "${OUTPUT_DIR}" \
        -t ${THREADS} \
        > "${log_file}" 2>&1; then
        
        echo "  ✓ Success"
        ((processed++))
        
        # Check if output files were created
        if [ -f "${sample_output}_report.tsv" ]; then
            clonotype_count=$(wc -l < "${sample_output}_report.tsv" | tr -d ' ')
            echo "  Clonotypes found: ${clonotype_count}"
        fi
    else
        echo "  ✗ Failed - see ${log_file} for details"
        ((failed++))
    fi
    
    echo ""
done

echo "=== Summary ==="
echo "Total samples: ${BAM_COUNT}"
echo "Processed: ${processed}"
echo "Skipped (already done): ${skipped}"
echo "Failed: ${failed}"
echo ""

# List output files
echo "=== Output Files ==="
ls -la "${OUTPUT_DIR}"/*.tsv 2>/dev/null | head -20 || echo "No TSV files found"

echo ""
echo "=== Next Steps ==="
echo "1. Run immunarch analysis:"
echo "   /Users/rishat/miniforge3/envs/multiomics/bin/Rscript \\"
echo "     ${BASE_DIR}/analysis/scripts/enhance_tcr_bcr_repertoire.R"
echo ""
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=== Complete ==="
