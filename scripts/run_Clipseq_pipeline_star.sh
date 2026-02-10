#!/usr/bin/env bash
set -euo pipefail

########################################
# CLIP/RIP-seq 自动化流程（更贴合CLIP特性）
# FastQC/MultiQC → fastp → (rRNA/tRNA预过滤) → STAR(genome) → 去重(可选) → featureCounts(可选)
#
# 修正点（相对你原脚本）：
# 1) 修复 SAMPLE 截取错误（_1/_2）与 DONE_FLAG 判断变量错误
# 2) fastp 参数对 CLIP 更温和（默认 q20 / len>=30；可按读长再调）
# 3) 真正使用你写好的 rRNA/tRNA STAR 索引进行预过滤（强烈建议）
# 4) STAR 多重比对策略与后续 unique bam 一致化（直接 Nmax=1）
# 5) 增加基础 mapping QC 输出与日志
# 6) 补上 featureCounts（你原脚本没跑）
########################################

# -------- 参数设置 --------
RAW_DIR="/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/Clip/CS2601026_人_clip/raw"
WORK_DIR="/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/new_all_res"

QC_DIR="$WORK_DIR/qc"
CLEAN_DIR="$WORK_DIR/clean_reads"
FILTER_DIR="$WORK_DIR/rrna_trna_filtered"
ALIGN_DIR="$WORK_DIR/alignment"
COUNT_DIR="$WORK_DIR/counts"

REF="/home/bio/biosoft/reference/GRCh38_UCSC/star_index"
GTF="/home/bio/biosoft/reference/GRCh38_UCSC/gtf/hg38.ncbiRefSeq.gtf"

RRNA_INDEX="/home/bio/biosoft/reference/GRCh38/rRNA_star_index"
TRNA_INDEX="/home/bio/biosoft/reference/GRCh38/tRNA_bt2_index/tRNA_dedup"

THREADS=32              # 96线程对IO/内存压力很大；建议按机器实际资源调整
FASTP_THREADS=16
STAR_THREADS=32
FEATURECOUNTS_THREADS=16

# 是否做 tRNA 预过滤（有些库 tRNA 很多，建议开；否则关也行）
DO_TRNA_FILTER=1

# 是否做去重：0=不去重；1=markdup（无UMI情形）；如果你是eCLIP带UMI，建议换成 umi-tools
DO_MARKDUP=1

# featureCounts：0=不跑；1=跑（CLIP 计数常用于QC/表达参考，不一定用于peak）
DO_FEATURECOUNTS=1

# fastp 更温和（更适合CLIP，避免过度丢读段）
FASTP_Q=20
FASTP_LEN=30

# STAR 基因组比对：如果你最终只用unique（-q255），这里直接限制multimap=1 更一致
STAR_MULTIMAP=1

# -------- 创建必要目录 --------
mkdir -p "$WORK_DIR" "$QC_DIR" "$CLEAN_DIR" "$FILTER_DIR" "$ALIGN_DIR" "$COUNT_DIR"

# -------- 工具检查 --------
for exe in fastqc multiqc fastp STAR samtools; do
  command -v "$exe" >/dev/null 2>&1 || { echo "❌ 缺少工具: $exe"; exit 1; }
done
if [[ "$DO_FEATURECOUNTS" -eq 1 ]]; then
  command -v featureCounts >/dev/null 2>&1 || { echo "❌ 缺少工具: featureCounts"; exit 1; }
fi

echo "✅ 参数概览:"
echo "RAW_DIR=$RAW_DIR"
echo "WORK_DIR=$WORK_DIR"
echo "THREADS: FastQC=$THREADS fastp=$FASTP_THREADS STAR=$STAR_THREADS FC=$FEATURECOUNTS_THREADS"
echo "fastp: q=$FASTP_Q len>=$FASTP_LEN"
echo "Filter: rRNA=ON tRNA=$DO_TRNA_FILTER"
echo "Markdup=$DO_MARKDUP featureCounts=$DO_FEATURECOUNTS"
echo

########################################
# Step 1: FastQC + MultiQC (raw)
########################################
echo "📊 Step 1: FastQC + MultiQC on raw fastq..."
QC_REPORT_RAW="$QC_DIR/raw/multiqc_report.html"
if [[ -f "$QC_REPORT_RAW" ]]; then
  echo "⏩ 已存在：$QC_REPORT_RAW"
else
  mkdir -p "$QC_DIR/raw"
  fastqc -o "$QC_DIR/raw" -t "$THREADS" "$RAW_DIR"/*fastq.gz
  export TZ="UTC"
  multiqc -o "$QC_DIR/raw" "$QC_DIR/raw"
  echo "✅ raw FastQC/MultiQC 完成：$QC_REPORT_RAW"
fi
echo

########################################
# Step 2: fastp 清洗（修正样本名截取）
########################################
echo "🧼 Step 2: fastp cleaning..."
shopt -s nullglob

for R1 in "$RAW_DIR"/*_1.fastq.gz; do
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  if [[ ! -f "$R2" ]]; then
    echo "⚠️ 找不到对应R2，跳过：$R1"
    continue
  fi

  BASE="$(basename "$R1")"
  SAMPLE="${BASE%_1.fastq.gz}"

  DONE_FLAG="$CLEAN_DIR/${SAMPLE}.done"
  if [[ -f "$DONE_FLAG" ]]; then
    echo "⏩ [$SAMPLE] fastp 已完成，跳过"
    continue
  fi

  echo "🧹 [$SAMPLE] fastp..."
  START=$(date +%s)

  fastp \
    -i "$R1" -I "$R2" \
    -o "$CLEAN_DIR/${SAMPLE}_clean_R1.fq.gz" \
    -O "$CLEAN_DIR/${SAMPLE}_clean_R2.fq.gz" \
    --detect_adapter_for_pe \
    -q "$FASTP_Q" \
    -u 30 \
    --length_required "$FASTP_LEN" \
    --cut_front --cut_tail \
    -h "$CLEAN_DIR/${SAMPLE}.fastp.html" \
    -j "$CLEAN_DIR/${SAMPLE}.fastp.json" \
    -w "$FASTP_THREADS" \
    > "$CLEAN_DIR/${SAMPLE}_fastp.log" 2>&1

  touch "$DONE_FLAG"
  END=$(date +%s)
  echo "✅ [$SAMPLE] fastp 完成：$((END - START)) 秒"
done

echo "🟢 fastp 完成，汇总 MultiQC..."
QC_REPORT_CLEAN="$CLEAN_DIR/multiqc_report.html"
if [[ -f "$QC_REPORT_CLEAN" ]]; then
  echo "⏩ 已存在：$QC_REPORT_CLEAN"
else
  export TZ="UTC"
  multiqc -o "$CLEAN_DIR" "$CLEAN_DIR"
  echo "✅ clean MultiQC 完成：$QC_REPORT_CLEAN"
fi
echo

########################################
# Step 3: rRNA（+可选tRNA）预过滤：把命中的reads丢掉，只保留unmapped进入主比对
# 说明：这里用 STAR 对 rRNA 索引跑一次，把未比对reads输出为 fastq，再进入 genome STAR。
########################################
echo "🧪 Step 3: rRNA/tRNA pre-filter with STAR..."
for R1 in "$CLEAN_DIR"/*_clean_R1.fq.gz; do
  R2="${R1/_clean_R1.fq.gz/_clean_R2.fq.gz}"
  BASE="$(basename "$R1")"
  SAMPLE="${BASE%_clean_R1.fq.gz}"

  OUT_R1="$FILTER_DIR/${SAMPLE}_nofilter_R1.fq.gz"
  OUT_R2="$FILTER_DIR/${SAMPLE}_nofilter_R2.fq.gz"

  DONE_FLAG="$FILTER_DIR/${SAMPLE}.done"
  if [[ -f "$DONE_FLAG" ]]; then
    echo "⏩ [$SAMPLE] 预过滤已完成，跳过"
    continue
  fi

  echo "🧪 [$SAMPLE] rRNA 预过滤..."
  TMPDIR="$FILTER_DIR/${SAMPLE}_tmp"
  mkdir -p "$TMPDIR"

  # rRNA
  STAR --runThreadN "$STAR_THREADS" \
    --genomeDir "$RRNA_INDEX" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand "pigz -dc" \
    --outFileNamePrefix "$TMPDIR/rrna_" \
    --outSAMtype None \
    --outReadsUnmapped Fastx \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMax 1 \
    > "$TMPDIR/rrna_filter.log" 2>&1

  # rRNA未比对reads（STAR输出为 Unmapped.out.mate1/2）
  RR1="$TMPDIR/rrna_Unmapped.out.mate1"
  RR2="$TMPDIR/rrna_Unmapped.out.mate2"

  # 可选 tRNA 再过滤
  if [[ "$DO_TRNA_FILTER" -eq 1 ]]; then
    echo "🧪 [$SAMPLE] tRNA 预过滤 (bowtie2)..."

    # bowtie2: 把“未命中 tRNA”的 paired reads 输出到指定文件
    # --un-conc-gz 会生成：prefix.1 / prefix.2（分别对应 mate1/mate2）
    # 这里用 OUT_R1/OUT_R2 作为最终进入 genome STAR 的 reads
    bowtie2 \
      -x "$TRNA_INDEX" \
      -1 "$RR1" -2 "$RR2" \
      --very-sensitive \
      -p "$STAR_THREADS" \
      --un-conc-gz "$TMPDIR/noTRNA.fq.gz" \
      -S /dev/null \
      > "$TMPDIR/trna_filter.log" 2>&1

    # bowtie2 输出命名规则：noTRNA_1.fq.gz / noTRNA_2.fq.gz
    TR1="$TMPDIR/noTRNA.fq.1.gz"
    TR2="$TMPDIR/noTRNA.fq.2.gz"

    # 复制/移动到标准输出名
    if [[ -f "$TR1" && -f "$TR2" ]]; then
      mv -f "$TR1" "$OUT_R1"
      mv -f "$TR2" "$OUT_R2"
    else
      echo "❌ [$SAMPLE] bowtie2 tRNA 过滤未产生预期输出：$TR1 / $TR2" | tee -a "$TMPDIR/trna_filter.log"
      exit 1
    fi

else
  # 不做 tRNA 过滤：直接把 rRNA 未命中 reads 作为输出
  pigz -c "$RR1" > "$OUT_R1"
  pigz -c "$RR2" > "$OUT_R2"
fi

  rm -rf "$TMPDIR"
  touch "$DONE_FLAG"
  echo "✅ [$SAMPLE] 预过滤完成：$OUT_R1 / $OUT_R2"
done
echo

########################################
# Step 4: STAR 基因组比对（更贴合CLIP：multimap策略一致化；输出基础QC）
########################################
echo "🧬 Step 4: STAR genome alignment..."
for R1 in "$FILTER_DIR"/*_nofilter_R1.fq.gz; do
  R2="${R1/_nofilter_R1.fq.gz/_nofilter_R2.fq.gz}"
  BASE="$(basename "$R1")"
  SAMPLE="${BASE%_nofilter_R1.fq.gz}"

  ALIGN_SAMPLE_DIR="$ALIGN_DIR/$SAMPLE"
  DONE_FLAG="$ALIGN_SAMPLE_DIR/${SAMPLE}.done"
  if [[ -f "$DONE_FLAG" ]]; then
    echo "⏩ [$SAMPLE] 比对已完成，跳过"
    continue
  fi

  mkdir -p "$ALIGN_SAMPLE_DIR"
  cd "$ALIGN_SAMPLE_DIR"

  ulimit -n 65535 || true

  echo "🧬 [$SAMPLE] STAR..."
  START=$(date +%s)

  STAR --runThreadN "$STAR_THREADS" \
    --genomeDir "$REF" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand "pigz -dc" \
    --outFileNamePrefix "${SAMPLE}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM NM MD \
    --outFilterType BySJout \
    --outFilterMultimapNmax "$STAR_MULTIMAP" \
    --outFilterMismatchNoverLmax 0.04 \
    --outFilterScoreMinOverLread 0.5 \
    --alignIntronMin 20 --alignIntronMax 500000 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --limitBAMsortRAM 30000000000 \
    > star.log 2>&1

  BAM="${SAMPLE}_Aligned.sortedByCoord.out.bam"

  # 如果你只要unique，这里输出 unique bam（q255）
  samtools view -@ "$STAR_THREADS" -b -q 255 "$BAM" > "${SAMPLE}.unique.bam"
  samtools index -@ "$STAR_THREADS" "${SAMPLE}.unique.bam"

  # 基础QC
  samtools flagstat -@ "$STAR_THREADS" "${SAMPLE}.unique.bam" > "${SAMPLE}.unique.flagstat.txt"
  samtools idxstats "${SAMPLE}.unique.bam" > "${SAMPLE}.unique.idxstats.txt"

  # 去重（可选；无UMI情况下用markdup，注意：需要先按name排序再fixmate再坐标排序）
  if [[ "$DO_MARKDUP" -eq 1 ]]; then
    echo "🧯 [$SAMPLE] markdup..."
    samtools sort -@ "$STAR_THREADS" -n -o "${SAMPLE}.nsort.bam" "${SAMPLE}.unique.bam"
    samtools fixmate -@ "$STAR_THREADS" -m "${SAMPLE}.nsort.bam" "${SAMPLE}.fixmate.bam"
    samtools sort -@ "$STAR_THREADS" -o "${SAMPLE}.fixmate.csort.bam" "${SAMPLE}.fixmate.bam"
    samtools markdup -@ "$STAR_THREADS" -r "${SAMPLE}.fixmate.csort.bam" "${SAMPLE}.unique.dedup.bam"
    samtools index -@ "$STAR_THREADS" "${SAMPLE}.unique.dedup.bam"
    samtools flagstat -@ "$STAR_THREADS" "${SAMPLE}.unique.dedup.bam" > "${SAMPLE}.unique.dedup.flagstat.txt"
    rm -f "${SAMPLE}.nsort.bam" "${SAMPLE}.fixmate.bam" "${SAMPLE}.fixmate.csort.bam"
  fi

  touch "$DONE_FLAG"
  END=$(date +%s)
  echo "✅ [$SAMPLE] 比对完成：$((END - START)) 秒"

  cd - >/dev/null
done
echo

########################################
# Step 5: featureCounts（可选）
# 说明：CLIP 的主分析通常是 peak calling；counts 更多用于QC/对照或表达参考。
# 建议用 dedup bam（如果做了去重），否则用 unique bam。
########################################
if [[ "$DO_FEATURECOUNTS" -eq 1 ]]; then
  echo "🧾 Step 5: featureCounts..."
  FC_OUT="$COUNT_DIR/featureCounts.gene.txt"
  if [[ -f "$FC_OUT" ]]; then
    echo "⏩ featureCounts 已存在：$FC_OUT"
  else
    BAM_LIST=()
    for d in "$ALIGN_DIR"/*; do
      [[ -d "$d" ]] || continue
      s="$(basename "$d")"
      if [[ "$DO_MARKDUP" -eq 1 && -f "$d/${s}.unique.dedup.bam" ]]; then
        BAM_LIST+=("$d/${s}.unique.dedup.bam")
      else
        BAM_LIST+=("$d/${s}.unique.bam")
      fi
    done

    # 这里按基因计数（-t exon -g gene_id 适用于常规GTF）
    featureCounts \
      -T "$FEATURECOUNTS_THREADS" \
      -a "$GTF" \
      -o "$FC_OUT" \
      -t exon -g gene_id \
      -p -B -C \
      "${BAM_LIST[@]}" \
      > "$COUNT_DIR/featureCounts.log" 2>&1

    echo "✅ featureCounts 完成：$FC_OUT"
  fi
fi

echo
echo "🎉 全流程完成！"
echo "1) raw QC:   $QC_DIR/raw/multiqc_report.html"
echo "2) clean QC: $CLEAN_DIR/multiqc_report.html"
echo "3) align:    $ALIGN_DIR/<sample>/"
echo "4) counts:   $COUNT_DIR/"
