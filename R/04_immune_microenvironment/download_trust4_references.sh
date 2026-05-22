#!/bin/bash
# ==============================================================================
# TRUST4 Reference Files Download Script
# ==============================================================================
# This script downloads the required reference files for TRUST4 analysis
# Run this script when network connectivity is available
# ==============================================================================

set -e

# Configuration
BASE_DIR="/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
REF_DIR="${BASE_DIR}/analysis/data/reference/trust4"

echo "=== TRUST4 Reference Files Download Script ==="
echo "Target directory: ${REF_DIR}"
echo ""

# Create directory
mkdir -p "${REF_DIR}"

# URLs for reference files
GITHUB_BASE="https://github.com/liulab-dfci/TRUST4/raw/master"
RAW_GITHUB="https://raw.githubusercontent.com/liulab-dfci/TRUST4/master"
JSDELIVR="https://cdn.jsdelivr.net/gh/liulab-dfci/TRUST4@master"
GHPROXY="https://ghp.ci/${GITHUB_BASE}"

# Files to download
FILES=("hg38_bcrtcr.fa" "human_IMGT+C.fa")

download_file() {
    local filename=$1
    local output_path="${REF_DIR}/${filename}"
    
    if [ -f "${output_path}" ] && [ -s "${output_path}" ]; then
        echo "✓ ${filename} already exists, skipping..."
        return 0
    fi
    
    echo "Downloading ${filename}..."
    
    # Try multiple sources
    local sources=(
        "${GITHUB_BASE}/${filename}"
        "${RAW_GITHUB}/${filename}"
        "${JSDELIVR}/${filename}"
        "${GHPROXY}/${filename}"
    )
    
    for url in "${sources[@]}"; do
        echo "  Trying: ${url}"
        if curl -L --connect-timeout 30 --max-time 300 -o "${output_path}" "${url}" 2>/dev/null; then
            # Check if file is valid (not an error page)
            if [ -s "${output_path}" ] && head -1 "${output_path}" | grep -q "^>"; then
                echo "  ✓ Success: Downloaded from ${url}"
                return 0
            else
                echo "  ✗ Invalid file content, trying next source..."
                rm -f "${output_path}"
            fi
        else
            echo "  ✗ Failed, trying next source..."
        fi
    done
    
    echo "  ✗ ERROR: Failed to download ${filename} from all sources"
    return 1
}

# Download each file
success_count=0
for file in "${FILES[@]}"; do
    if download_file "${file}"; then
        ((success_count++))
    fi
done

echo ""
echo "=== Download Summary ==="
echo "Successfully downloaded: ${success_count}/${#FILES[@]} files"

if [ ${success_count} -eq ${#FILES[@]} ]; then
    echo ""
    echo "✓ All reference files downloaded successfully!"
    echo "  - hg38_bcrtcr.fa: Human T/B cell receptor reference (hg38)"
    echo "  - human_IMGT+C.fa: IMGT reference with constant regions"
    echo ""
    echo "You can now run TRUST4 analysis with:"
    echo "  bash ${BASE_DIR}/analysis/scripts/run_trust4_analysis.sh"
else
    echo ""
    echo "✗ Some files failed to download."
    echo "Please check your network connection and try again."
    echo "Alternatively, manually download from:"
    echo "  ${GITHUB_BASE}/hg38_bcrtcr.fa"
    echo "  ${GITHUB_BASE}/human_IMGT+C.fa"
fi

echo ""
echo "=== Complete ==="
