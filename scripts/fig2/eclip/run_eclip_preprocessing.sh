#!/usr/bin/env bash
set -euo pipefail

########################################
# Fig. 2b eCLIP preprocessing and peak calling
# (STAR genome alignment, keep multimappers)
# 1) FastQC/MultiQC (raw)
# 2) fastp (company-like: -u 50, explicit adapters, len>=30)
# 3) rRNA filter (STAR rRNA index) -> keep unmapped reads
# 4) STAR genome alignment (ALLOW multimapping, NOT unique-only)
# 5) Build MACS2 BAM:
#    - primary-only (+ optional proper pair)
#    - PLUS TLEN upper bound filter to avoid absurd fragment size in BAMPE
# 6) Optional Picard MarkDuplicates (REMOVE_DUPLICATES=true) on TLEN-filtered BAM
# 7) MACS2 callpeak (BAMPE, --nomodel, -p 0.01) + tighter peak width control
########################################

# ---------------- USER CONFIG ----------------
PROJECT_ROOT="${XAF1_PROJECT_ROOT:-$(pwd)}"
REFERENCE_ROOT="${XAF1_REFERENCE_ROOT:?Set XAF1_REFERENCE_ROOT to the directory containing the genome indexes}"
RAW_DIR="${XAF1_ECLIP_RAW_DIR:-${PROJECT_ROOT}/data/eclip/raw}"
WORK_DIR="${XAF1_ECLIP_WORK_DIR:-${PROJECT_ROOT}/results/fig2/eclip_preprocessing}"

QC_DIR="$WORK_DIR/qc"
CLEAN_DIR="$WORK_DIR/clean_reads"
FILTER_DIR="$WORK_DIR/rrna_filtered"
ALIGN_DIR="$WORK_DIR/alignment"
MACS2_DIR="$WORK_DIR/macs2"
TMP_DIR="$WORK_DIR/tmp"   # <-- NEW: shared tmp dir for STAR

# STAR indexes
REF="${XAF1_STAR_INDEX:-${REFERENCE_ROOT}/GRCh38_UCSC/star_index}"
RRNA_INDEX="${XAF1_RRNA_STAR_INDEX:-${REFERENCE_ROOT}/GRCh38/rRNA_star_index}"

# Threads
THREADS="${XAF1_THREADS:-16}"
FASTQC_THREADS="${XAF1_FASTQC_THREADS:-${THREADS}}"
FASTP_THREADS="${XAF1_FASTP_THREADS:-${THREADS}}"
STAR_THREADS="${XAF1_STAR_THREADS:-${THREADS}}"
SAMTOOLS_THREADS="${XAF1_SAMTOOLS_THREADS:-${THREADS}}"

# fastp (company-like)
FASTP_Q=20
FASTP_U=50
FASTP_LEN=30
ADAPTER_R1="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATC"

# rRNA filter
DO_RRNA_FILTER=1
RRNA_STAR_THREADS=64

# STAR multimapping policy (KEEP multimappers)
STAR_MULTIMAP=30

# Build BAM for MACS2:
MACS2_REQUIRE_PROPER_PAIR=1

# NEW: TLEN upper bound for BAMPE sanity
# If your library is typical CLIP/eCLIP, TLEN rarely > 1000-2000.
# Start with 2000. If too strict/too many reads lost, relax to 3000-5000.
MACS2_TLEN_MAX=1500

# Optional dedup (company uses Picard; you can turn it on/off here)
DO_PICARD_DEDUP=1
PICARD_JAR=""   # leave empty if `picard` command exists; else set /path/to/picard.jar

# MACS2
GENOME_SIZE="hs"
MACS2_PVALUE="0.01"
MACS2_QVALUE="0.05"
MAKE_TRACKS=1   # adds -B --SPMR

# Peak width control (tighten to avoid overly wide peaks)
MACS2_MAX_GAP=15
MACS2_MIN_LENGTH=30

# Broad peaks tend to be wide; default OFF
DO_BROAD=0
BROAD_CUTOFF="0.1"

# bigWig generation (requires UCSC tools)
MAKE_BW=0
CHROM_SIZES="${XAF1_CHROM_SIZES:-${REFERENCE_ROOT}/GRCh38_UCSC/fasta/hg38.chrom.sizes}"
BDG_THREADS="${XAF1_BDG_THREADS:-${THREADS}}"

# Per-sample bigWig from BAM (independent of IP/Input)
MAKE_SAMPLE_BW=1
BW_NORMALIZE="CPM"   # CPM=每百万reads；也可用 RPKM / RPGC(需effectiveGenomeSize)
BW_BINSIZE=10        # 10bp bins 常用；想更平滑可用 25
BW_THREADS="${XAF1_BW_THREADS:-${THREADS}}"

ulimit -n 65535 || true

# ---------------- END CONFIG ----------------

log(){ echo -e "[$(date '+%F %T')] $*"; }

mkdir -p "$QC_DIR"/{raw,clean} "$CLEAN_DIR" "$FILTER_DIR" "$ALIGN_DIR" \
         "$MACS2_DIR"/{logs,narrow,broad,tmp,tracks} "$TMP_DIR"

# ---------------- Tool checks ----------------
for exe in fastqc multiqc fastp STAR samtools macs2 java awk; do
  command -v "$exe" >/dev/null 2>&1 || { echo "❌ Missing tool: $exe"; exit 1; }
done

run_picard() {
  if command -v picard >/dev/null 2>&1; then
    picard "$@"
  else
    if [[ -z "$PICARD_JAR" ]]; then
      echo "❌ picard not found in PATH and PICARD_JAR is empty."
      exit 1
    fi
    java -jar "$PICARD_JAR" "$@"
  fi
}

if [[ "${MAKE_BW:-0}" -eq 1 ]]; then
  command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "❌ Missing tool: bedGraphToBigWig"; exit 1; }
fi

mkdir -p "$MACS2_DIR/tracks/samples"

if [[ "${MAKE_SAMPLE_BW:-0}" -eq 1 ]]; then
  command -v bamCoverage >/dev/null 2>&1 || { echo "❌ Missing tool: bamCoverage (deeptools)"; exit 1; }
fi

log "✅ Config summary:"
log "RAW_DIR=$RAW_DIR"
log "WORK_DIR=$WORK_DIR"
log "fastp: q=$FASTP_Q u=$FASTP_U len>=$FASTP_LEN"
log "rRNA filter: $DO_RRNA_FILTER"
log "STAR: threads=$STAR_THREADS multimapNmax=$STAR_MULTIMAP (keep multimappers)"
log "MACS2 BAM: primary-only; proper-pair=$MACS2_REQUIRE_PROPER_PAIR; TLEN_MAX=$MACS2_TLEN_MAX"
log "Picard dedup: $DO_PICARD_DEDUP"
log "MACS2: -p $MACS2_PVALUE --nomodel --max-gap $MACS2_MAX_GAP --min-length $MACS2_MIN_LENGTH"
echo

########################################
# Step 1: FastQC + MultiQC (raw)
########################################
log "������ Step 1: FastQC + MultiQC on raw..."
RAW_MQC="$QC_DIR/raw/multiqc_report.html"
if [[ -f "$RAW_MQC" ]]; then
  log "⏩ Exists: $RAW_MQC"
else
  fastqc -o "$QC_DIR/raw" -t "$FASTQC_THREADS" "$RAW_DIR"/*fastq.gz
  multiqc -o "$QC_DIR/raw" "$QC_DIR/raw"
  log "✅ Raw QC done: $RAW_MQC"
fi
echo

########################################
# Step 2: fastp cleaning
########################################
log "������ Step 2: fastp cleaning..."
shopt -s nullglob
for R1 in "$RAW_DIR"/*_1.fastq.gz; do
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  [[ -f "$R2" ]] || { log "⚠️ Missing R2 for $R1 (skip)"; continue; }

  BASE="$(basename "$R1")"
  SAMPLE="${BASE%_1.fastq.gz}"
  DONE="$CLEAN_DIR/${SAMPLE}.done"
  [[ -f "$DONE" ]] && { log "⏩ [$SAMPLE] fastp done"; continue; }

  log "������ [$SAMPLE] fastp..."
  fastp \
    -i "$R1" -I "$R2" \
    -o "$CLEAN_DIR/${SAMPLE}_clean_R1.fq.gz" \
    -O "$CLEAN_DIR/${SAMPLE}_clean_R2.fq.gz" \
    --adapter_sequence "$ADAPTER_R1" \
    --adapter_sequence_r2 "$ADAPTER_R2" \
    -q "$FASTP_Q" \
    -u "$FASTP_U" \
    --length_required "$FASTP_LEN" \
    --cut_front --cut_tail \
    -w "$FASTP_THREADS" \
    -h "$CLEAN_DIR/${SAMPLE}.fastp.html" \
    -j "$CLEAN_DIR/${SAMPLE}.fastp.json" \
    > "$CLEAN_DIR/${SAMPLE}.fastp.log" 2>&1

  touch "$DONE"
  log "✅ [$SAMPLE] fastp done"
done

CLEAN_MQC="$QC_DIR/clean/multiqc_report.html"
if [[ -f "$CLEAN_MQC" ]]; then
  log "⏩ Exists: $CLEAN_MQC"
else
  multiqc -o "$QC_DIR/clean" "$CLEAN_DIR"
  log "✅ Clean QC done: $CLEAN_MQC"
fi
echo

########################################
# Step 3: rRNA prefilter (STAR) -> keep unmapped reads
########################################
log "������ Step 3: rRNA prefilter..."
for R1 in "$CLEAN_DIR"/*_clean_R1.fq.gz; do
  R2="${R1/_clean_R1.fq.gz/_clean_R2.fq.gz}"
  BASE="$(basename "$R1")"
  SAMPLE="${BASE%_clean_R1.fq.gz}"

  OUT_R1="$FILTER_DIR/${SAMPLE}_rrnaFree_R1.fq.gz"
  OUT_R2="$FILTER_DIR/${SAMPLE}_rrnaFree_R2.fq.gz"
  DONE="$FILTER_DIR/${SAMPLE}.done"
  [[ -f "$DONE" ]] && { log "⏩ [$SAMPLE] rRNA filter done"; continue; }

  if [[ "$DO_RRNA_FILTER" -eq 0 ]]; then
    ln -sf "$R1" "$OUT_R1"
    ln -sf "$R2" "$OUT_R2"
    touch "$DONE"
    log "ℹ️ [$SAMPLE] rRNA filter OFF (symlink clean reads)"
    continue
  fi

  TMP="$FILTER_DIR/${SAMPLE}_tmp"
  mkdir -p "$TMP"

  log "������ [$SAMPLE] STAR rRNA filter..."
  STAR --runThreadN "$RRNA_STAR_THREADS" \
    --genomeDir "$RRNA_INDEX" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand "pigz -dc" \
    --outFileNamePrefix "$TMP/rrna_" \
    --outSAMtype None \
    --outReadsUnmapped Fastx \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMax 1 \
    --outTmpDir "$TMP_DIR/${SAMPLE}.rrna_STARtmp" \
    > "$TMP/rrna_filter.log" 2>&1

  pigz -c "$TMP/rrna_Unmapped.out.mate1" > "$OUT_R1"
  pigz -c "$TMP/rrna_Unmapped.out.mate2" > "$OUT_R2"

  rm -rf "$TMP" "$TMP_DIR/${SAMPLE}.rrna_STARtmp" 2>/dev/null || true
  touch "$DONE"
  log "✅ [$SAMPLE] rRNA-free reads ready"
done
echo

########################################
# Step 4: STAR genome alignment (keep multimappers)
########################################
log "������ Step 4: STAR genome alignment (multimappers kept)..."
for R1 in "$FILTER_DIR"/*_rrnaFree_R1.fq.gz; do
  R2="${R1/_rrnaFree_R1.fq.gz/_rrnaFree_R2.fq.gz}"
  BASE="$(basename "$R1")"
  SAMPLE="${BASE%_rrnaFree_R1.fq.gz}"

  SDIR="$ALIGN_DIR/$SAMPLE"
  DONE="$SDIR/${SAMPLE}.star.done"
  [[ -f "$DONE" ]] && { log "⏩ [$SAMPLE] STAR done"; continue; }

  mkdir -p "$SDIR"
  pushd "$SDIR" >/dev/null

  log "������ [$SAMPLE] STAR..."
  STAR --runThreadN "$STAR_THREADS" \
    --genomeDir "$REF" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand "pigz -dc" \
    --outFileNamePrefix "${SAMPLE}_" \
    --outTmpDir "$TMP_DIR/${SAMPLE}.genome_STARtmp" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM NM MD \
    --outFilterMultimapNmax "$STAR_MULTIMAP" \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 500000 \
    --limitBAMsortRAM 30000000000 \
    > star.log 2>&1

  BAM="${SAMPLE}_Aligned.sortedByCoord.out.bam"
  samtools index -@ "$SAMTOOLS_THREADS" "$BAM"
  samtools flagstat -@ "$SAMTOOLS_THREADS" "$BAM" > "${SAMPLE}.star.flagstat.txt"

  rm -rf "$TMP_DIR/${SAMPLE}.genome_STARtmp" 2>/dev/null || true

  touch "$DONE"
  popd >/dev/null
  log "✅ [$SAMPLE] STAR alignment done"
done
echo

########################################
# Step 5: Build MACS2 input BAM (primary-only + proper pair + TLEN cap)
########################################
log "������ Step 5: Build MACS2 BAM (primary-only, keep multimappers, TLEN-capped)..."
for d in "$ALIGN_DIR"/*; do
  [[ -d "$d" ]] || continue
  SAMPLE="$(basename "$d")"
  IN_BAM="$d/${SAMPLE}_Aligned.sortedByCoord.out.bam"
  [[ -f "$IN_BAM" ]] || continue

  # keep your original output too (for debugging)
  OUT_PRIMARY="$d/${SAMPLE}.macs2.primary.bam"

  # NEW: final BAM for MACS2 (after TLEN filter)
  OUT_PP_TLEN="$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.bam"
  DONE="$d/${SAMPLE}.macs2bam.done"
  [[ -f "$DONE" ]] && { log "⏩ [$SAMPLE] MACS2 BAM exists"; continue; }

  log "������ [$SAMPLE] filtering primary alignments..."
  if [[ "$MACS2_REQUIRE_PROPER_PAIR" -eq 1 ]]; then
    samtools view -@ "$SAMTOOLS_THREADS" -b -f 2 -F 256 -F 2048 "$IN_BAM" \
      | samtools sort -@ "$SAMTOOLS_THREADS" -o "$OUT_PRIMARY" -
  else
    samtools view -@ "$SAMTOOLS_THREADS" -b 10 -F 256 -F 2048 "$IN_BAM" \
      | samtools sort -@ "$SAMTOOLS_THREADS" -o "$OUT_PRIMARY" -
  fi
  samtools index -@ "$SAMTOOLS_THREADS" "$OUT_PRIMARY"
  samtools flagstat -@ "$SAMTOOLS_THREADS" "$OUT_PRIMARY" > "$d/${SAMPLE}.macs2.primary.flagstat.txt"

  log "������ [$SAMPLE] applying TLEN cap (<=${MACS2_TLEN_MAX}) for BAMPE sanity..."
  # Keep header; for each alignment, compute abs(TLEN) from field 9, keep 0<TLEN<=MAX
  # Note: in proper PE, TLEN should be >0 or <0 depending on read; abs() handles both.
  samtools view -h "$OUT_PRIMARY" \
  | awk -v MAX="$MACS2_TLEN_MAX" 'BEGIN{OFS="\t"} /^@/{print; next} {t=$9; if(t<0) t=-t; if(t>0 && t<=MAX) print}' \
  | samtools view -@ "$SAMTOOLS_THREADS" -b - \
  | samtools sort -@ "$SAMTOOLS_THREADS" -o "$OUT_PP_TLEN" -

  samtools index -@ "$SAMTOOLS_THREADS" "$OUT_PP_TLEN"
  samtools flagstat -@ "$SAMTOOLS_THREADS" "$OUT_PP_TLEN" > "$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.flagstat.txt"

  # Optional: quick TLEN quantiles to a text file (helps choose MAX)
  samtools view -@ "$SAMTOOLS_THREADS" -f 2 -F 256 -F 2048 "$OUT_PP_TLEN" \
  | awk '{t=$9; if(t<0) t=-t; if(t>0) print t}' \
  | sort -n \
  | awk '{
      a[NR]=$1
    } END {
      if(NR>0){
        print "n="NR, "p50="a[int(NR*0.50)], "p90="a[int(NR*0.90)], "p99="a[int(NR*0.99)], "max="a[NR]
      } else {
        print "n=0"
      }
    }' > "$d/${SAMPLE}.tlen_stats.after_cap.txt" || true

  touch "$DONE"
  log "✅ [$SAMPLE] MACS2 input BAM: $OUT_PP_TLEN"
done
echo

########################################
# Step 6: Optional Picard MarkDuplicates (remove duplicates) on TLEN-capped BAM
########################################
if [[ "$DO_PICARD_DEDUP" -eq 1 ]]; then
  log "������ Step 6: Picard MarkDuplicates (REMOVE_DUPLICATES=true)..."
  for d in "$ALIGN_DIR"/*; do
    [[ -d "$d" ]] || continue
    SAMPLE="$(basename "$d")"

    IN_BAM="$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.bam"
    [[ -f "$IN_BAM" ]] || continue

    OUT_BAM="$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.dedup.bam"
    METRICS="$d/${SAMPLE}.picard.dup_metrics.txt"
    DONE="$d/${SAMPLE}.picard.done"
    [[ -f "$DONE" ]] && { log "⏩ [$SAMPLE] Picard done"; continue; }

    log "������ [$SAMPLE] Picard MarkDuplicates..."
    run_picard MarkDuplicates \
      I="$IN_BAM" \
      O="$OUT_BAM" \
      M="$METRICS" \
      REMOVE_DUPLICATES=true \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      > "$d/${SAMPLE}.picard.log" 2>&1

    samtools index -@ "$SAMTOOLS_THREADS" "$OUT_BAM"
    samtools flagstat -@ "$SAMTOOLS_THREADS" "$OUT_BAM" > "$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.dedup.flagstat.txt"

    touch "$DONE"
    log "✅ [$SAMPLE] dedup BAM: $OUT_BAM"
  done
  echo
fi

########################################
# Step 7: MACS2 callpeak (BAMPE, --nomodel, -p 0.01, width control)
########################################
log "������ Step 7: MACS2 callpeak..."

sample_key(){ echo "${1#*-}"; }

mapfile -t IP_SAMPLES < <(ls -1 "$ALIGN_DIR" 2>/dev/null | grep '^IP-' || true)
mapfile -t IN_SAMPLES < <(ls -1 "$ALIGN_DIR" 2>/dev/null | grep '^Input-' || true)

if [[ "${#IP_SAMPLES[@]}" -eq 0 ]]; then
  echo "❌ No IP-* folders found under: $ALIGN_DIR"; exit 1
fi
if [[ "${#IN_SAMPLES[@]}" -eq 0 ]]; then
  echo "❌ No Input-* folders found under: $ALIGN_DIR"; exit 1
fi

declare -A INPUT_BY_KEY
for s in "${IN_SAMPLES[@]}"; do
  INPUT_BY_KEY["$(sample_key "$s")"]="$s"
done

pick_call_bam() {
  local dir="$1"
  local sample="$2"
  # Prefer TLEN-capped dedup if enabled and exists
  if [[ "$DO_PICARD_DEDUP" -eq 1 && -f "$dir/${sample}.macs2.pp.tlen${MACS2_TLEN_MAX}.dedup.bam" ]]; then
    echo "$dir/${sample}.macs2.pp.tlen${MACS2_TLEN_MAX}.dedup.bam"
  elif [[ -f "$dir/${sample}.macs2.pp.tlen${MACS2_TLEN_MAX}.bam" ]]; then
    echo "$dir/${sample}.macs2.pp.tlen${MACS2_TLEN_MAX}.bam"
  else
    echo ""
  fi
}

for ip in "${IP_SAMPLES[@]}"; do
  key="$(sample_key "$ip")"
  input="${INPUT_BY_KEY[$key]:-}"
  [[ -n "$input" ]] || { log "⚠️ No matching Input for $ip (key=$key), skip"; continue; }

  IP_DIR="$ALIGN_DIR/$ip"
  IN_DIR="$ALIGN_DIR/$input"

  IP_BAM="$(pick_call_bam "$IP_DIR" "$ip")"
  IN_BAM="$(pick_call_bam "$IN_DIR" "$input")"

  [[ -n "$IP_BAM" ]] || { log "❌ Missing IP MACS2 BAM for $ip"; continue; }
  [[ -n "$IN_BAM" ]] || { log "❌ Missing Input MACS2 BAM for $input"; continue; }

  PREFIX="${ip}_vs_${input}"
  NDIR="$MACS2_DIR/narrow/$PREFIX"
  mkdir -p "$NDIR"

  log "������ Pair: IP=$ip  Input=$input"
  log "   IP_BAM=$IP_BAM"
  log "   IN_BAM=$IN_BAM"

  if [[ ! -f "$NDIR/.done" ]]; then
    log "������ MACS2 narrow (control width): $PREFIX"
    macs2 callpeak \
      -t "$IP_BAM" \
      -c "$IN_BAM" \
      -f BAMPE \
      -g "$GENOME_SIZE" \
      -p "$MACS2_PVALUE" \
      --nomodel \
      --call-summits \
      --keep-dup all \
      --max-gap "$MACS2_MAX_GAP" \
      --min-length "$MACS2_MIN_LENGTH" \
      $( [[ "$MAKE_TRACKS" -eq 1 ]] && echo "-B --SPMR" ) \
      -n "$PREFIX" \
      --outdir "$NDIR" \
      > "$MACS2_DIR/logs/${PREFIX}.narrow.log" 2>&1
    touch "$NDIR/.done"
    log "✅ Narrow done: $NDIR"
  else
    log "⏩ Narrow exists: $PREFIX"
  fi

  if [[ "$DO_BROAD" -eq 1 ]]; then
    BDIR="$MACS2_DIR/broad/$PREFIX"
    mkdir -p "$BDIR"
    if [[ ! -f "$BDIR/.done" ]]; then
      log "������ MACS2 broad: $PREFIX"
      macs2 callpeak \
        -t "$IP_BAM" \
        -c "$IN_BAM" \
        -f BAMPE \
        -g "$GENOME_SIZE" \
        -p "$MACS2_PVALUE" \
        --nomodel \
        --broad --broad-cutoff "$BROAD_CUTOFF" \
        --keep-dup all \
        $( [[ "$MAKE_TRACKS" -eq 1 ]] && echo "-B --SPMR" ) \
        -n "$PREFIX" \
        --outdir "$BDIR" \
        > "$MACS2_DIR/logs/${PREFIX}.broad.log" 2>&1
      touch "$BDIR/.done"
      log "✅ Broad done: $BDIR"
    else
      log "⏩ Broad exists: $PREFIX"
    fi
  fi
done

########################################
# Step 8: bedGraph -> bigWig (optional)
# Requires: -B --SPMR in MACS2 (MAKE_TRACKS=1) and CHROM_SIZES + bedGraphToBigWig
########################################
if [[ "${MAKE_BW:-0}" -eq 1 ]]; then
  log "������ Step 8: Convert MACS2 bedGraph -> bigWig..."

  if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "❌ CHROM_SIZES not found: $CHROM_SIZES"
    echo "   Please set CHROM_SIZES to a valid hg38.chrom.sizes"
    exit 1
  fi

  # Helper: sanitize and sort bdg for bedGraphToBigWig
  bdg_to_bw() {
    local bdg="$1"
    local bw="$2"

    # MACS2 sometimes outputs unsorted or has out-of-bound lines; sanitize:
    # 1) remove header/track lines
    # 2) keep only lines with start<end
    # 3) sort by chr, start
    awk 'BEGIN{OFS="\t"} $1!~/^track/ && $1!~/^browser/ && NF>=4 {if($2<$3) print $1,$2,$3,$4}' "$bdg" \
      | sort -k1,1 -k2,2n --parallel="${BDG_THREADS:-4}" \
      > "${bdg}.sorted"

    bedGraphToBigWig "${bdg}.sorted" "$CHROM_SIZES" "$bw"
    rm -f "${bdg}.sorted"
  }

  # Iterate all MACS2 narrow outputs and convert treat/control bdg and log2FE if present
  for ndir in "$MACS2_DIR"/narrow/*; do
    [[ -d "$ndir" ]] || continue
    prefix="$(basename "$ndir")"

    # MACS2 outputs names like: <name>_treat_pileup.bdg / <name>_control_lambda.bdg
    treat="$ndir/${prefix}_treat_pileup.bdg"
    ctrl="$ndir/${prefix}_control_lambda.bdg"

    out_treat_bw="$MACS2_DIR/tracks/${prefix}.treat.SPMR.bw"
    out_ctrl_bw="$MACS2_DIR/tracks/${prefix}.control.lambda.bw"

    if [[ -f "$treat" ]]; then
      if [[ ! -f "$out_treat_bw" ]]; then
        log "  ������ treat bdg -> bw: $prefix"
        bdg_to_bw "$treat" "$out_treat_bw"
      else
        log "  ⏩ treat bw exists: $out_treat_bw"
      fi
    else
      log "  ⚠️ treat bdg missing (did MACS2 run with -B --SPMR?): $treat"
    fi

    if [[ -f "$ctrl" ]]; then
      if [[ ! -f "$out_ctrl_bw" ]]; then
        log "  ������ control bdg -> bw: $prefix"
        bdg_to_bw "$ctrl" "$out_ctrl_bw"
      else
        log "  ⏩ control bw exists: $out_ctrl_bw"
      fi
    else
      log "  ⚠️ control bdg missing: $ctrl"
    fi

    # Optional: make log2 fold-enrichment track (requires both bdg)
    # MACS2 provides bdgcmp to generate logFE
    if command -v macs2 >/dev/null 2>&1 && [[ -f "$treat" && -f "$ctrl" ]]; then
      logfe_bdg="$MACS2_DIR/tracks/${prefix}.log2FE.bdg"
      logfe_bw="$MACS2_DIR/tracks/${prefix}.log2FE.bw"
      if [[ ! -f "$logfe_bw" ]]; then
        log "  ������ bdgcmp log2FE: $prefix"
        macs2 bdgcmp -t "$treat" -c "$ctrl" -m logFE -o "$logfe_bdg" \
          > "$MACS2_DIR/logs/${prefix}.bdgcmp.log" 2>&1 || true
        if [[ -f "$logfe_bdg" ]]; then
          bdg_to_bw "$logfe_bdg" "$logfe_bw"
          rm -f "$logfe_bdg"
        fi
      else
        log "  ⏩ log2FE bw exists: $logfe_bw"
      fi
    fi
  done

  log "✅ bigWig tracks saved to: $MACS2_DIR/tracks/"
fi

########################################
# Step 9: Per-sample bigWig from BAM (IP and Input separately)
# Uses deepTools bamCoverage. This is NOT IP-vs-Input; it's per-sample coverage.
########################################
if [[ "${MAKE_SAMPLE_BW:-0}" -eq 1 ]]; then
  log "������ Step 9: Generate per-sample bigWig from BAM (independent tracks)..."

  for d in "$ALIGN_DIR"/*; do
    [[ -d "$d" ]] || continue
    SAMPLE="$(basename "$d")"

    # Prefer TLEN-capped dedup BAM (matches MACS2 input cleanliness)
    BAM=""
    if [[ "$DO_PICARD_DEDUP" -eq 1 && -f "$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.dedup.bam" ]]; then
      BAM="$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.dedup.bam"
    elif [[ -f "$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.bam" ]]; then
      BAM="$d/${SAMPLE}.macs2.pp.tlen${MACS2_TLEN_MAX}.bam"
    elif [[ -f "$d/${SAMPLE}_Aligned.sortedByCoord.out.bam" ]]; then
      # fallback: raw STAR BAM
      BAM="$d/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    else
      continue
    fi

    # Ensure BAM is indexed
    [[ -f "${BAM}.bai" ]] || samtools index -@ "$SAMTOOLS_THREADS" "$BAM"

    OUTBW="$MACS2_DIR/tracks/samples/${SAMPLE}.${BW_NORMALIZE}.bs${BW_BINSIZE}.bw"
    if [[ -f "$OUTBW" ]]; then
      log "  ⏩ bw exists: $OUTBW"
      continue
    fi

    log "  ������ bamCoverage -> bw: $SAMPLE"
    bamCoverage \
      -b "$BAM" \
      -o "$OUTBW" \
      --binSize "$BW_BINSIZE" \
      --normalizeUsing "$BW_NORMALIZE" \
      -p "$BW_THREADS" \
      > "$MACS2_DIR/logs/${SAMPLE}.bamCoverage.log" 2>&1

    log "  ✅ saved: $OUTBW"
  done

  log "✅ Per-sample bigWigs saved to: $MACS2_DIR/tracks/samples/"
fi


echo
log "������ All done!"
log "Raw QC:   $QC_DIR/raw/multiqc_report.html"
log "Clean QC: $QC_DIR/clean/multiqc_report.html"
log "Align:    $ALIGN_DIR/<sample>/"
log "MACS2:    $MACS2_DIR/"
