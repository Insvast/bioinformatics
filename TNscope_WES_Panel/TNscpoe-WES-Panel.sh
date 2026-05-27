#!/bin/bash
set -xeuo pipefail

if [ "$#" -lt 11 ]; then
  echo "Usage: $0 PAIR_ID WORKDIR TUMOR_SM TUMOR_R1 TUMOR_R2 NORMAL_SM NORMAL_R1 NORMAL_R2 FASTA DBSNP INTERVAL_FILE [TNSCOPE_FILTER]" >&2
  exit 1
fi

PAIR_ID=$1
WORKDIR=$2
TUMOR_SM=$3
TUMOR_FASTQ_1=$4
TUMOR_FASTQ_2=$5
NORMAL_SM=$6
NORMAL_FASTQ_1=$7
NORMAL_FASTQ_2=$8
FASTA=$9
KNOWN_DBSNP=${10}
INTERVAL_FILE=${11}
TNSCOPE_FILTER=${12:-tnscope_filter.py}
PLATFORM=${PLATFORM:-auto}

TUMOR_RGID="rg_${TUMOR_SM}"
NORMAL_RGID="rg_${NORMAL_SM}"

unset http_proxy
unset https_proxy

source /miniconda3/etc/profile.d/conda.sh
conda activate sentieon_cli_env

export SENTIEON_INSTALL_DIR=/APP/u22/x86_com/sentieon/202503.03/sentieon-genomics-202503.03
export SENTIEON_LICENSE="需要指定"
export PATH=$SENTIEON_INSTALL_DIR/bin:$PATH
export LD_LIBRARY_PATH=$SENTIEON_INSTALL_DIR/lib:${LD_LIBRARY_PATH:-}
export LD_PRELOAD=/miniconda3/envs/anasys/lib/libjemalloc.so.2

NT=$(nproc)
LOGFILE=${PAIR_ID}.run.log

for f in "$TUMOR_FASTQ_1" "$TUMOR_FASTQ_2" "$NORMAL_FASTQ_1" "$NORMAL_FASTQ_2" "$FASTA" "$KNOWN_DBSNP"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: required input not found or empty: $f" >&2
    exit 1
  fi
done

if [ -n "$INTERVAL_FILE" ] && [ "$INTERVAL_FILE" != "NA" ] && [ ! -s "$INTERVAL_FILE" ]; then
  echo "ERROR: interval BED not found or empty: $INTERVAL_FILE" >&2
  exit 1
fi

if [ ! -s "$TNSCOPE_FILTER" ]; then
  echo "ERROR: TNscope filter script not found or empty: $TNSCOPE_FILTER" >&2
  exit 1
fi

if [ "$PLATFORM" = "auto" ]; then
  read_fields=$(zcat "$TUMOR_FASTQ_1" 2>/dev/null | head -n 1 | awk '{print NF}' || true)
  if [ -z "$read_fields" ]; then
    echo "ERROR: failed to read first FASTQ record from $TUMOR_FASTQ_1" >&2
    exit 1
  fi
  if [ "$read_fields" -eq 1 ]; then
    PLATFORM=DNBSEQ
  else
    PLATFORM=ILLUMINA
  fi
fi

mkdir -p "$WORKDIR"
cd "$WORKDIR"
exec >"$LOGFILE" 2>&1

echo "PairID: $PAIR_ID"
echo "Tumor: $TUMOR_SM $TUMOR_FASTQ_1 $TUMOR_FASTQ_2"
echo "Normal: $NORMAL_SM $NORMAL_FASTQ_1 $NORMAL_FASTQ_2"
echo "FASTA: $FASTA"
echo "DBSNP: $KNOWN_DBSNP"
echo "INTERVAL_FILE: $INTERVAL_FILE"
echo "PLATFORM: $PLATFORM"
echo "THREADS: $NT"
echo "Start time: $(date "+%Y-%m-%d %H:%M:%S")"

interval_arg=()
if [ -n "$INTERVAL_FILE" ] && [ "$INTERVAL_FILE" != "NA" ]; then
  interval_arg=(--interval "$INTERVAL_FILE")
fi

(
  "$SENTIEON_INSTALL_DIR/bin/sentieon" bwa mem -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PLATFORM" \
    -t "$NT" -K 10000000 "$FASTA" "$TUMOR_FASTQ_1" "$TUMOR_FASTQ_2" || \
    { echo "BWA tumor error"; exit 1; }
) | "$SENTIEON_INSTALL_DIR/bin/sentieon" util sort -o tumor_sorted.bam -t "$NT" --sam2bam -i - || \
  { echo "Tumor alignment failed"; exit 1; }

(
  "$SENTIEON_INSTALL_DIR/bin/sentieon" bwa mem -R "@RG\tID:$NORMAL_RGID\tSM:$NORMAL_SM\tPL:$PLATFORM" \
    -t "$NT" -K 10000000 "$FASTA" "$NORMAL_FASTQ_1" "$NORMAL_FASTQ_2" || \
    { echo "BWA normal error"; exit 1; }
) | "$SENTIEON_INSTALL_DIR/bin/sentieon" util sort -o normal_sorted.bam -t "$NT" --sam2bam -i - || \
  { echo "Normal alignment failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -r "$FASTA" -t "$NT" -i tumor_sorted.bam \
  --algo MeanQualityByCycle tumor_mq_metrics.txt \
  --algo QualDistribution tumor_qd_metrics.txt --algo GCBias \
  --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat \
  --adapter_seq '' tumor_aln_metrics.txt \
  --algo InsertSizeMetricAlgo tumor_is_metrics.txt \
  --algo CoverageMetrics --omit_base_output tumor_coverage_metrics || \
  { echo "Tumor metrics failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -r "$FASTA" -t "$NT" -i normal_sorted.bam \
  --algo MeanQualityByCycle normal_mq_metrics.txt \
  --algo QualDistribution normal_qd_metrics.txt --algo GCBias \
  --summary normal_gc_summary.txt normal_gc_metrics.txt --algo AlignmentStat \
  --adapter_seq '' normal_aln_metrics.txt \
  --algo InsertSizeMetricAlgo normal_is_metrics.txt \
  --algo CoverageMetrics --omit_base_output normal_coverage_metrics || \
  { echo "Normal metrics failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -t "$NT" -i tumor_sorted.bam --algo LocusCollector \
  --fun score_info tumor_score.txt || { echo "Tumor LocusCollector failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -t "$NT" -i tumor_sorted.bam --algo Dedup \
  --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam || \
  { echo "Tumor Dedup failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -t "$NT" -i normal_sorted.bam --algo LocusCollector \
  --fun score_info normal_score.txt || { echo "Normal LocusCollector failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -t "$NT" -i normal_sorted.bam --algo Dedup \
  --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam || \
  { echo "Normal Dedup failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -r "$FASTA" -t "$NT" \
  "${interval_arg[@]}" \
  -i tumor_deduped.bam -i normal_deduped.bam \
  --algo TNscope \
  --tumor_sample "$TUMOR_SM" --normal_sample "$NORMAL_SM" \
  --dbsnp "$KNOWN_DBSNP" \
  --min_tumor_allele_frac 0.009 \
  --max_fisher_pv_active 0.05 \
  --max_normal_alt_cnt 10 \
  --max_normal_alt_frac 0.01 \
  --max_normal_alt_qsum 250 \
  --sv_mask_ext 10 \
  --prune_factor 0 \
  output_tnscope.pre_filter.vcf.gz || \
  { echo "TNscope failed"; exit 1; }

"$SENTIEON_INSTALL_DIR/bin/sentieon" pyexec "$TNSCOPE_FILTER" \
  -v output_tnscope.pre_filter.vcf.gz \
  --tumor_sample "$TUMOR_SM" \
  --normal_sample "$NORMAL_SM" \
  -x tissue_panel --min_tumor_af 0.0095 --min_depth 200 \
  output_tnscope.filter.vcf.gz || \
  { echo "TNscope filter failed"; exit 1; }

echo "End time: $(date "+%Y-%m-%d %H:%M:%S")"
