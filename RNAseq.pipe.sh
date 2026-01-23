#!/bin/bash
#SBATCH --job-name=WL-ERRXXXXXXXX
#SBATCH --partition=deimos
#SBATCH --cpus-per-task=64

# ==============================================================================
# Sentieon RNA-seq 变异检测流程脚本
# 基于 Sentieon Release 202503 文档 - Chapter 6
# ==============================================================================


#!/bin/bash
echo $0 \$SAMPLEID \$WORKDIR \$FASTQ_1 \$FASTQ_2 \$FASTA \$STAR_INDEX 
set -euxo pipefail

export SAMPLEID=$1
export WORKDIR=$2
export FASTQ_1=$3
export FASTQ_2=$4
export REF_FASTA=$5
export STAR_INDEX=$6

LOGFILE=$SAMPLEID.run.log
export root=/mnt/Bio/softwore/sentieon-genomics-202503.02
export SENTIEON_INSTALL_DIR="/mnt/Bio/softwore/sentieon-genomics-202503.02"
SENTIEON_EXEC="$SENTIEON_INSTALL_DIR/bin/sentieon"
export NT=$(nproc)  # 64


# 样本信息 (用于 Read Group @RG)
GROUP_NAME=rg_${SAMPLEID}
SAMPLE_NAME=${SAMPLEID}
PLATFORM="ILLUMINA"
READ_LENGTH_MINUS_1=149
# 输出目录
WORKDIR=${WORKDIR}/result/${SAMPLEID}
mkdir -p $WORKDIR
cd $WORKDIR
exec >>$LOGFILE 2>&1
export clean1=clean.$(basename $FASTQ_1)
export clean2=clean.$(basename $FASTQ_2)

timer(){
    start_time=$(date +%s)
    start_date==$(date -d @$start_time "+%Y-%m-%d %H:%M:%S")
    eval $2 && touch $3
    end_time=$(date +%s)
    start_date==$(date -d @$end_time "+%Y-%m-%d %H:%M:%S")
    cost_time=$[ $end_time-$start_time ]
    echo -e "TIMER: $1\t$(($cost_time/60)) min $(($cost_time%60)) s"
}


raw2clean(){
    cmd="fastp -w 16 -i $FASTQ_1 -I $FASTQ_2 -o $clean1 -O $clean2 -j $SAMPLEID.qc.json -h $SAMPLEID.qc.html \
    --detect_adapter_for_pe \
    --length_required 30 -q 15 -u 40 -n 5 \
    --trim_front1 0 --trim_tail1 0 --trim_front2 0 --trim_tail2 0"
    timer raw2clean "$cmd" qc.ok
}

alignment(){
    cmd="$SENTIEON_EXEC STAR \
        --runThreadN $NT \
        --genomeDir $STAR_INDEX \
        --readFilesIn $clean1 $clean2 \
        --readFilesCommand "zcat" \
        --outStd BAM_Unsorted \
        --outSAMtype BAM Unsorted \
        --outBAMcompression 0 \
        --twopass1readsN -1 \
        --outSAMattrRGline ID:${GROUP_NAME} SM:${SAMPLE_NAME} PL:${PLATFORM} \
        --twopassMode Basic \
        --sjdbOverhang $READ_LENGTH_MINUS_1 \
        | $SENTIEON_EXEC util sort -r $REF_FASTA -o ${SAMPLE_NAME}.sorted.bam -t $NT -i -"
    timer alignment "$cmd" align.ok
}

dedup(){
    cmd="$SENTIEON_EXEC driver -t $NT -i ${SAMPLE_NAME}.sorted.bam \
        --algo LocusCollector \
        --rna \
        --fun score_info score.gz"
    timer LocusCollector "$cmd" LocusCollector.ok

    cmd="$SENTIEON_EXEC driver -t $NT -i ${SAMPLE_NAME}.sorted.bam \
        --algo Dedup \
        --score_info score.gz \
        --metrics dedup_metrics.txt \
        ${SAMPLE_NAME}.deduped.bam"
    timer Dedup "$cmd" dedup.ok
}

split(){
    cmd="$SENTIEON_EXEC driver -t $NT -r $REF_FASTA -i ${SAMPLE_NAME}.deduped.bam \
    --algo RNASplitReadsAtJunction \
    --reassign_mapq 255:60 \
    ${SAMPLE_NAME}.split.bam"
    timer RNASplitReadsAtJunction "$cmd" split.ok
}

DNAscope(){
    cmd="$SENTIEON_EXEC driver -t $NT -r $REF_FASTA -i ${SAMPLE_NAME}.split.bam \
    --algo DNAscope --trim_soft_clip  \
    --call_conf 20 --emit_conf 20 ${SAMPLE_NAME}.rna.vcf.gz"
    timer DNAscope "$cmd" dnascope.ok
}


metrics(){
    # 统计质控指标
    cmd="$SENTIEON_EXEC driver -t $NT -r $REF_FASTA -i ${SAMPLE_NAME}.sorted.bam \
        --algo MeanQualityByCycle $SAMPLEID.mq_metrics.txt \
        --algo QualDistribution $SAMPLEID.qd_metrics.txt   \
        --algo GCBias --summary $SAMPLEID.gc_summary.txt $SAMPLEID.gc_metrics.txt \
        --algo AlignmentStat --adapter_seq ''  $SAMPLEID.aln_metrics.txt  \
        --algo InsertSizeMetricAlgo $SAMPLEID.is_metrics.txt"
    timer metrics "$cmd" metrics.ok
}


[ -e qc.ok ]||raw2clean
[ -e qc.ok ]&&[ ! -e align.ok ]&&alignment
[ -e dedup.ok ]||dedup
[ -e split.ok ]||split
[ -e dnascope.ok ]||DNAscope
[ -e metrics.ok ]||metrics

echo ">>> Pipeline Finished Successfully!"
echo ">>> 结束时间：$(date "+%Y-%m-%d %H:%M:%S")"
