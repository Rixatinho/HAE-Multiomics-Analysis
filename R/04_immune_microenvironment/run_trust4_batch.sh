#!/bin/bash
# TRUST4 batch processing with space-safe path handling

cd "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"

BAM_DIR="/Volumes/阿热-实验数据备份/工作/1.肝包虫/1.泡型包虫/1.2024-11-泡球蚴/1.转录组学测序/02.Bam"
OUT_DIR="analysis/results/enhancement_tcr_bcr/trust4_output"
REF_FA="analysis/data/reference/trust4/hg38_bcrtcr.fa"
IMGT_FA="analysis/data/reference/trust4/human_IMGT+C.fa"
LOG_DIR="analysis/results/enhancement_tcr_bcr/logs"
TRUST4="/Users/rishat/miniforge3/envs/multiomics/bin/run-trust4"

mkdir -p "$OUT_DIR" "$LOG_DIR"

processed=0
skipped=0
failed=0
total=26

echo "=== TRUST4 Batch Processing ==="
echo "Start: $(date)"
echo ""

for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    
    # Skip if already processed
    if [ -f "$OUT_DIR/${sample}_report.tsv" ]; then
        echo "[SKIP] $sample - already processed"
        ((skipped++))
        continue
    fi
    
    echo "[RUN] $sample ($((processed + skipped + failed + 1))/$total)"
    
    # Run TRUST4 with relative paths
    if "$TRUST4" \
        -b "$bam" \
        -f "$REF_FA" \
        --ref "$IMGT_FA" \
        -o "$sample" \
        --od "$OUT_DIR" \
        -t 4 \
        > "$LOG_DIR/${sample}.log" 2>&1; then
        
        echo "  ✓ Success"
        if [ -f "$OUT_DIR/${sample}_report.tsv" ]; then
            count=$(wc -l < "$OUT_DIR/${sample}_report.tsv" | tr -d ' ')
            echo "  Clonotypes: $count"
        fi
        ((processed++))
    else
        echo "  ✗ Failed"
        ((failed++))
    fi
done

echo ""
echo "=== Summary ==="
echo "Processed: $processed"
echo "Skipped: $skipped"
echo "Failed: $failed"
echo "End: $(date)"
