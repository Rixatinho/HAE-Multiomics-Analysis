#!/bin/bash
# 包装脚本：运行 nc_production_pipeline.R 并捕获输出

cd "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"

# 激活 conda 环境
eval "$(conda shell.zsh hook)"
conda activate multiomics

# 运行 pipeline
echo "=== Starting nc_production_pipeline.R at $(date) ===" > pipeline_output.log
Rscript analysis/scripts/nc_production_pipeline.R >> pipeline_output.log 2>&1
EXIT_CODE=$?

echo "" >> pipeline_output.log
echo "=== Pipeline finished at $(date) with exit code $EXIT_CODE ===" >> pipeline_output.log

# 统计生成的文件
echo "" >> pipeline_output.log
echo "=== Output files in 图片/ subdirectories ===" >> pipeline_output.log
find "图片" -maxdepth 2 -name "*.pdf" -path "*/Figure_*/*" -o -name "*.pdf" -path "*/SuppFig_*/*" 2>/dev/null | wc -l >> pipeline_output.log
