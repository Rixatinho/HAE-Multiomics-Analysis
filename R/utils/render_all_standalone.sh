#!/bin/bash
# 批量渲染所有30个standalone脚本
# 使用方法: bash render_all_standalone.sh

set -e  # 不要因为单个失败就停止

# 工作目录
WORKDIR="/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学"
cd "$WORKDIR"

# 日志文件
LOG_FILE="$WORKDIR/render_standalone_log.txt"
SUMMARY_FILE="$WORKDIR/render_standalone_summary.txt"

# 清空日志
> "$LOG_FILE"
> "$SUMMARY_FILE"

echo "========================================" | tee -a "$LOG_FILE"
echo "开始批量渲染 $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

# 计数器
SUCCESS=0
FAILED=0
TOTAL=0

# 函数：运行单个脚本
run_script() {
    local script_path="$1"
    local script_name=$(basename "$script_path")
    local dir_name=$(dirname "$script_path")
    
    TOTAL=$((TOTAL + 1))
    echo "" | tee -a "$LOG_FILE"
    echo "[$TOTAL/30] 正在运行: $script_name" | tee -a "$LOG_FILE"
    echo "开始时间: $(date '+%H:%M:%S')" | tee -a "$LOG_FILE"
    
    START_TIME=$(date +%s)
    
    # 运行脚本 (使用miniforge3绝对路径)
    if /Users/rishat/miniforge3/envs/multiomics/bin/Rscript "$script_path" >> "$LOG_FILE" 2>&1; then
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))
        echo "✓ 成功: $script_name (耗时: ${DURATION}s)" | tee -a "$LOG_FILE"
        echo "SUCCESS: $script_name (${DURATION}s)" >> "$SUMMARY_FILE"
        SUCCESS=$((SUCCESS + 1))
    else
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))
        echo "✗ 失败: $script_name (耗时: ${DURATION}s)" | tee -a "$LOG_FILE"
        echo "FAILED: $script_name (${DURATION}s)" >> "$SUMMARY_FILE"
        FAILED=$((FAILED + 1))
    fi
}

# 运行8张主图
for i in $(seq 1 8); do
    SCRIPT="图片/Figure_$i/Figure_${i}_standalone.R"
    if [ -f "$SCRIPT" ]; then
        run_script "$SCRIPT"
    else
        echo "警告: 脚本不存在 - $SCRIPT" | tee -a "$LOG_FILE"
        echo "NOT_FOUND: $SCRIPT" >> "$SUMMARY_FILE"
        FAILED=$((FAILED + 1))
        TOTAL=$((TOTAL + 1))
    fi
done

# 运行22张补充图
for i in $(seq 1 22); do
    NUM=$(printf "%02d" $i)
    SCRIPT="图片/SuppFig_$NUM/SuppFig_${NUM}_standalone.R"
    if [ -f "$SCRIPT" ]; then
        run_script "$SCRIPT"
    else
        echo "警告: 脚本不存在 - $SCRIPT" | tee -a "$LOG_FILE"
        echo "NOT_FOUND: $SCRIPT" >> "$SUMMARY_FILE"
        FAILED=$((FAILED + 1))
        TOTAL=$((TOTAL + 1))
    fi
done

# 汇总
echo "" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "渲染完成 $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "总计: $TOTAL 个脚本" | tee -a "$LOG_FILE"
echo "成功: $SUCCESS 个" | tee -a "$LOG_FILE"
echo "失败: $FAILED 个" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

# 写入汇总文件末尾
echo "" >> "$SUMMARY_FILE"
echo "TOTAL: $TOTAL" >> "$SUMMARY_FILE"
echo "SUCCESS_COUNT: $SUCCESS" >> "$SUMMARY_FILE"
echo "FAILED_COUNT: $FAILED" >> "$SUMMARY_FILE"
