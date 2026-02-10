#!/usr/bin/env bash
set -euo pipefail

#############################################
# MACS2 peak calling for CLIP (paired-end BAMPE)
# - Auto-pair IP-* with Input-* by suffix after the first dash
# - Prefer *.unique.dedup.bam, fallback to *.unique.bam
# - Run both broad and narrow+summits
# - Optional peaklets (summits ±75bp)
# - Optional bdg->bw (requires bedGraphToBigWig + chrom.sizes)
#
# Author: ChatGPT + Shiyou Wang
#############################################

# --------- USER CONFIG ---------
ALIGN_ROOT="/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/new_all_res/alignment"
OUT_ROOT="/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/new_all_res/macs2"

THREADS=16
GENOME_SIZE="hs"   # hs for human; or use effective genome size number
Q_NARROW="0.05"
Q_BROAD="0.1"
BROAD_CUTOFF="0.1"

# For peaklets & bw generation:
MAKE_PEAKLETS=1
PEAKLET_HALF_WIDTH=75  # summit ±75bp => 150bp peaklets

MAKE_BW=1
CHROM_SIZES="/home/bio/biosoft/reference/GRCh38_UCSC/fasta/hg38.chrom.sizes"  # <-- set if MAKE_BW=1

# --------- TOOL CHECK ---------
for exe in macs2 samtools bedtools; do
  command -v "$exe" >/dev/null 2>&1 || { echo "❌ Missing tool: $exe"; exit 1; }
done
if [[ "$MAKE_BW" -eq 1 ]]; then
  command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "❌ Missing tool: bedGraphToBigWig"; exit 1; }
fi

mkdir -p "$OUT_ROOT"/{broad,narrow,peaklets,logs,tracks,tmp}

# --------- HELPERS ---------
pick_bam() {
  local sample_dir="$1"
  local sample="$2"
  local bam1="${sample_dir}/${sample}.unique.dedup.bam"
  local bam2="${sample_dir}/${sample}.unique.bam"
  if [[ -f "$bam1" ]]; then
    echo "$bam1"
  elif [[ -f "$bam2" ]]; then
    echo "$bam2"
  else
    echo ""
  fi
}

# Pairing rule:
#   IP-1-EV71-1  -> key="1-EV71-1"
#   Input-1-EV71-1 -> same key
sample_key() {
  local s="$1"
  echo "${s#*-}"   # strip up to first '-'
}

log() {
  echo -e "[$(date '+%F %T')] $*"
}

# --------- DISCOVER SAMPLES ---------
log "🔎 Scanning alignment folders: $ALIGN_ROOT"

mapfile -t IP_SAMPLES < <(ls -1 "$ALIGN_ROOT" | grep '^IP-' || true)
mapfile -t IN_SAMPLES < <(ls -1 "$ALIGN_ROOT" | grep '^Input-' || true)

if [[ "${#IP_SAMPLES[@]}" -eq 0 ]]; then
  echo "❌ No IP-* folders found under: $ALIGN_ROOT"
  exit 1
fi
if [[ "${#IN_SAMPLES[@]}" -eq 0 ]]; then
  echo "❌ No Input-* folders found under: $ALIGN_ROOT"
  exit 1
fi

# Build associative map for Input by key
declare -A INPUT_BY_KEY
for s in "${IN_SAMPLES[@]}"; do
  k="$(sample_key "$s")"
  INPUT_BY_KEY["$k"]="$s"
done

# --------- MAIN LOOP ---------
for ip in "${IP_SAMPLES[@]}"; do
  k="$(sample_key "$ip")"
  input="${INPUT_BY_KEY[$k]:-}"

  if [[ -z "$input" ]]; then
    log "⚠️  No matching Input for $ip (key=$k). Skipping."
    continue
  fi

  IP_DIR="$ALIGN_ROOT/$ip"
  IN_DIR="$ALIGN_ROOT/$input"

  IP_BAM="$(pick_bam "$IP_DIR" "$ip")"
  IN_BAM="$(pick_bam "$IN_DIR" "$input")"

  if [[ -z "$IP_BAM" ]]; then
    log "❌ Missing BAM for IP: $ip (expected *.unique(.dedup).bam). Skipping."
    continue
  fi
  if [[ -z "$IN_BAM" ]]; then
    log "❌ Missing BAM for Input: $input (expected *.unique(.dedup).bam). Skipping."
    continue
  fi

  log "🧬 Pair: IP=$ip  Input=$input"
  log "   IP_BAM=$IP_BAM"
  log "   IN_BAM=$IN_BAM"

  # Quick sanity: BAM indexed?
  if [[ ! -f "${IP_BAM}.bai" ]]; then samtools index -@ "$THREADS" "$IP_BAM"; fi
  if [[ ! -f "${IN_BAM}.bai" ]]; then samtools index -@ "$THREADS" "$IN_BAM"; fi

  PREFIX="${ip}_vs_${input}"

  # ---------- Broad calling ----------
  BROAD_DIR="$OUT_ROOT/broad/$PREFIX"
  mkdir -p "$BROAD_DIR"
  DONE_BROAD="$BROAD_DIR/.done"
  if [[ -f "$DONE_BROAD" ]]; then
    log "⏩ Broad already done: $PREFIX"
  else
    log "🚀 MACS2 broad calling: $PREFIX"
    macs2 callpeak \
      -t "$IP_BAM" \
      -c "$IN_BAM" \
      -f BAMPE \
      -g "$GENOME_SIZE" \
      --broad --broad-cutoff "$BROAD_CUTOFF" \
      -q "$Q_BROAD" \
      --keep-dup all \
      -n "$PREFIX" \
      --outdir "$BROAD_DIR" \
      > "$OUT_ROOT/logs/${PREFIX}.broad.log" 2>&1
    touch "$DONE_BROAD"
    log "✅ Broad done: $BROAD_DIR"
  fi

  # ---------- Narrow + summits ----------
  NARROW_DIR="$OUT_ROOT/narrow/$PREFIX"
  mkdir -p "$NARROW_DIR"
  DONE_NARROW="$NARROW_DIR/.done"
  if [[ -f "$DONE_NARROW" ]]; then
    log "⏩ Narrow already done: $PREFIX"
  else
    log "🚀 MACS2 narrow+summits: $PREFIX"
    macs2 callpeak \
      -t "$IP_BAM" \
      -c "$IN_BAM" \
      -f BAMPE \
      -g "$GENOME_SIZE" \
      -q "$Q_NARROW" \
      --call-summits \
      --keep-dup all \
      -n "$PREFIX" \
      --outdir "$NARROW_DIR" \
      > "$OUT_ROOT/logs/${PREFIX}.narrow.log" 2>&1
    touch "$DONE_NARROW"
    log "✅ Narrow done: $NARROW_DIR"
  fi

  # ---------- Peaklets (summits ±W) ----------
  if [[ "$MAKE_PEAKLETS" -eq 1 ]]; then
    SUMMITS="$NARROW_DIR/${PREFIX}_summits.bed"
    if [[ -f "$SUMMITS" ]]; then
      PEAKLET_OUT="$OUT_ROOT/peaklets/${PREFIX}.peaklets_${PEAKLET_HALF_WIDTH}bp.bed"
      if [[ -f "$PEAKLET_OUT" ]]; then
        log "⏩ Peaklets exist: $PEAKLET_OUT"
      else
        log "🧩 Building peaklets (±${PEAKLET_HALF_WIDTH}bp): $PREFIX"
        # Requires chrom sizes for bedtools slop. If you don't have CHROM_SIZES here,
        # you can generate it from your reference or from the STAR genome.
        if [[ ! -f "$CHROM_SIZES" ]]; then
          echo "❌ CHROM_SIZES not found: $CHROM_SIZES"
          echo "   Set CHROM_SIZES to hg38.chrom.sizes or disable MAKE_PEAKLETS."
          exit 1
        fi
        bedtools slop -i "$SUMMITS" -g "$CHROM_SIZES" -b "$PEAKLET_HALF_WIDTH" \
          | sort -k1,1 -k2,2n > "$PEAKLET_OUT"
        log "✅ Peaklets saved: $PEAKLET_OUT"
      fi
    else
      log "⚠️ Summits not found (narrow may have failed?): $SUMMITS"
    fi
  fi

  # ---------- Optional: tracks (bdg/bw) ----------
  if [[ "$MAKE_BW" -eq 1 ]]; then
    if [[ ! -f "$CHROM_SIZES" ]]; then
      echo "❌ CHROM_SIZES not found: $CHROM_SIZES"
      exit 1
    fi

    # Use MACS2 bdgcmp for log2(IP/Input) tracks from pileup bedGraphs if present.
    # MACS2 callpeak may output *_treat_pileup.bdg and *_control_lambda.bdg when -B is used.
    # If you want tracks, you should add -B --SPMR to callpeak commands above.
    log "ℹ️  To generate bw tracks, enable -B --SPMR in macs2 callpeak and re-run."
  fi

done

log "🎉 All MACS2 peak calling done."
log "Results:"
log "  Broad:   $OUT_ROOT/broad"
log "  Narrow:  $OUT_ROOT/narrow"
log "  Peaklets:$OUT_ROOT/peaklets"
