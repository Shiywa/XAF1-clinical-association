#!/bin/bash

########################################
# 自动化 RIP-seq 流程：FastQC → fastp → STAR + featureCounts
# 作者: ChatGPT + Shiyou Wang
########################################

# -------- 参数设置 --------
RAW_DIR="/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/Clip/CS2601026_人_clip/raw"        
WORK_DIR="/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res"
CLEAN_DIR="$WORK_DIR/clean_reads"
ALIGN_DIR="$WORK_DIR/alignment"
QC_DIR="$WORK_DIR/qc"

REF="/home/bio/biosoft/reference/GRCh38_UCSC/star_index"
GTF="/home/bio/biosoft/reference/GRCh38_UCSC/gtf/hg38.ncbiRefSeq.gtf"
RRNA_INDEX=/home/bio/biosoft/reference/GRCh38/rRNA_star_index        # 你已构建好的 rRNA STAR 索引目录
tRNA_INDEX=/home/bio/biosoft/reference/GRCh38/tRNA_star_index        # 你已构建好的 tRNA STAR 索引目录

THREADS=96
featurecount_thread=64
FASTP_THREADS=96

use_R=/home/bio/micromamba/envs/R440/bin/Rscript

# -------- 创建必要目录 --------
mkdir -p $WORK_DIR $QC_DIR $CLEAN_DIR $ALIGN_DIR

# -------- Step 1: FastQC + MultiQC --------
echo "📊 Step 1: Running FastQC..."
QC_REPORT="$QC_DIR/raw/multiqc_report.html"

if [ -f "$QC_REPORT" ]; then
    echo "✅ QA报告已存在，跳过FastQC与MultiQC分析：$QC_REPORT"
else
    mkdir -p $QC_DIR/raw
    fastqc -o $QC_DIR/raw -t $THREADS $RAW_DIR/*fastq.gz
    export TZ="UTC"
    multiqc -o $QC_DIR/raw $QC_DIR/raw
    echo "✅ FastQC与MultiQC分析完成！"
fi

# -------- Step 2: fastp 清洗 --------
echo "🧼 Step 2: Cleaning reads using fastp..."
for F in $RAW_DIR/*_1.fastq.gz; do
    START=$(date +%s)

    R=${F/_1/_2}
    BASE=${F##*/}
    SAMPLE=${BASE%_2.fastq.gz}

    DONE_FLAG="$CLEAN_DIR/${SAMPLE}.done"
    if [ -f "$DONE_FLAG" ]; then
        echo "⏩ [$SAMPLE] 已完成清洗，跳过"
        continue
    fi

    echo "🧹 Cleaning $SAMPLE ..."
    if fastp       -i $F -I $R       -o $CLEAN_DIR/${SAMPLE}_clean_R1.fq.gz       -O $CLEAN_DIR/${SAMPLE}_clean_R2.fq.gz       --detect_adapter_for_pe       -q 25       -u 20       -e 20       -r       -W 5       -M 30       --length_required 50       -h $CLEAN_DIR/${SAMPLE}.html       -j $CLEAN_DIR/${SAMPLE}.json       -w $FASTP_THREADS > $CLEAN_DIR/${SAMPLE}_fastp.log 2>&1; then
        touch "$DONE_FLAG"
        END=$(date +%s)
        echo "✅ [$SAMPLE] Fastp 清洗完成！耗时 $((END - START)) 秒"
    else
        echo "❌ [$SAMPLE] Fastp 失败" | tee "$CLEAN_DIR/${SAMPLE}.error"
        continue
    fi
done

echo "🟢 所有样本 fastp 清洗完成，开始整合 QC 报告..."

QC_REPORT2="$CLEAN_DIR/multiqc_report.html"

if [ -f "$QC_REPORT2" ]; then
    echo "✅ fastp 清洗报告已存在，跳过fastp 清洗和MultiQC分析：$QC_REPORT2"
else
    start_time=$(date +%s)
    export TZ="UTC"
    multiqc $CLEAN_DIR -o $CLEAN_DIR
    end_time=$(date +%s)
    echo "✅ fastp 清洗与MultiQC分析完成！ (用时 $((end_time - start_time)) 秒)"
fi

# -------- Step 3:  --------
echo "🧬 Step 3: Aligning and quantifying with STAR + featureCounts..."
for F in $CLEAN_DIR/*_clean_R1.fq.gz; do
    START=$(date +%s)

    BASE=${F##*/}
    SAMPLE=${BASE%_clean_R1.fq.gz}

    ALIGN_SAMPLE_DIR=$ALIGN_DIR/$SAMPLE
    DONE_FLAG="$ALIGN_SAMPLE_DIR/${SAMPLE}.done"

    if [ -f "$DONE_FLAG1" ]; then
        echo "⏩ [$SAMPLE] 已完成比对，跳过"
        continue
    fi

    echo "🧬 Processing $SAMPLE ..."
    mkdir -p $ALIGN_SAMPLE_DIR
    cd $ALIGN_SAMPLE_DIR

    ulimit -n 65535

    # Step 1: STAR 比对
    if STAR --runThreadN $THREADS \
            --genomeDir $REF \
            --readFilesIn $CLEAN_DIR/${SAMPLE}_clean_R1.fq.gz $CLEAN_DIR/${SAMPLE}_clean_R2.fq.gz \
            --readFilesCommand "pigz -dc" \
	    --outFileNamePrefix ${SAMPLE}_ \
	    --outSAMtype BAM SortedByCoordinate \
	    --outSAMunmapped Within \
	    --outSAMattributes NH HI AS nM NM MD \
	    --outFilterType BySJout \
	    --outFilterMultimapNmax 10 \
	    --outFilterMismatchNoverLmax 0.04 \
	    --outFilterScoreMinOverLread 0.5 \
	    --alignIntronMin 20 --alignIntronMax 500000 \
	    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
	    --sjdbScore 1 > star.log 2>&1; then
        echo "✅ [$SAMPLE] STAR比对完成！"

	samtools view -b -q 255 ${SAMPLE}_Aligned.sortedByCoord.out.bam > ${SAMPLE}_unique.bam

	# 生成标志文件
        touch "$DONE_FLAG"
        END=$(date +%s)
        echo "✅ [$SAMPLE] 比对与定量完成！耗时 $((END - START)) 秒"
    else
        echo "❌ [$SAMPLE] 比对失败" | tee "${SAMPLE}.error"
    fi

    cd - > /dev/null
done

