#!/bin/bash
echo $0 \$SAMPLEID \$WORKDIR \$FASTQ_1 \$FASTQ_2 \$FASTA \$DATATYPE \$KEEP_CLEAN \$KEEP_BAM \$PL \$BED \$PCR_FREE
# set -euxo pipefail

SAMPLEID=$1
WORKDIR=$2
FASTQ_1=$3
FASTQ_2=$4
FASTA=$5
TYPE=${6:-"raw"}
KEEP_CLEAN=${7:-keep}
KEEP_BAM=${8:-keep}
PL=${9:-2}
BED=${10:-"none"}
PCR_FREE=${11:-"false"}
LOGFILE=$SAMPLEID.run.log
SENTIEON_INSTALL_DIR=/sentieon/202503.03/sentieon-genomics-202503.03
module load bcftools/1.21-gcc-11.4.0 samtools/1.20-gcc-11.4.0
SENTIEON_EXEC="${SENTIEON_INSTALL_DIR}/bin/sentieon"
LD_PRELOAD=/miniconda3/envs/anasys/lib/libjemalloc.so.2
# 需要配置
SENTIEON_LICENSE=IP:8091

# 主程序开始
# 判断测序平台并使用对应模型参数
if [ -e $LOGFILE ];then
    export PLATFORM=$(awk '/Platform/{print $NF;exit}' $LOGFILE)
else
    export count=$(zcat $FASTQ_1|head -n 1|awk '{print NF}')
    if [ $count -eq 1 ];then
        export PLATFORM="DNBSEQ"
    elif [ $count == 2 ];then
        export PLATFORM="ILLUMINA"
    elif [ $count == 3 ];then
        export PLATFORM="ILLUMINA"
    else
        echo "Unrecognized platform"
        export PLATFORM="ILLUMINA"
    fi
fi
echo $PLATFORM

# 文件后缀
FAI=$FASTA.fai

export SUFFIX=$(awk 'BEGIN{max=0}{if($2>max)max=$2}END{if(max>536870911){print "cram"}else{print "bam"}}' $FAI)

export TMPDIR=$WORKDIR

# 分析线程数
export THREADS=$(nproc)

# 创建分析目录
[ -e $WORKDIR ]||mkdir -p $WORKDIR
cd $WORKDIR
exec >>$LOGFILE 2>&1
echo "SampleID:" $SAMPLEID
echo "DataType:" $TYPE
echo "Platform:" $PLATFORM
echo "THREADS:" $THREADS
echo "KEEP_CLEAN:" $KEEP_CLEAN
echo "KEEP_BAM:" $KEEP_BAM
echo "TMPDIR:" $TMPDIR

export bwt_max_mem=200G
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
    cmd="fastp -w 16 -i $FASTQ_1 -I $FASTQ_2 -o $clean1 -O $clean2 -j $SAMPLEID.qc.json -h $SAMPLEID.qc.html"
    timer raw2clean "$cmd" qc.ok
}



dnascope(){
    tag="@RG\tID:rg_$SAMPLEID\tSM:$SAMPLEID\tPL:$PLATFORM"
    cmd="sentieon-cli dnascope -r $FASTA --r1-fastq $clean1 --r2-fastq $clean2 --readgroups \"$tag\" -m $ML_MODEL $bed_arg --interval_padding 100 -t $THREADS -g --duplicate-marking markdup --assay $type_arg --consensus $pcr_arg $suffix_arg $SAMPLEID.vcf.gz && md5sum $SAMPLEID.vcf.gz > $SAMPLEID.vcf.gz.md5"
    timer dnascope "$cmd" dnascope.ok
}


saferm(){
    for i in $*
    do
        [ -e $i ]&&rm -rf $i||echo $i
    done
}

if [ $TYPE == "raw" ];then
    export clean1=clean.$(basename $FASTQ_1)
    export clean2=clean.$(basename $FASTQ_2)
    [ -e qc.ok ]||raw2clean
else
    export clean1=$FASTQ_1
    export clean2=$FASTQ_2
    [ -e qc.ok ]||touch qc.ok
fi

if [ $PLATFORM == "DNBSEQ" ];then
    export ML_MODEL=/resource/DNAscopeMGIWGS2.0.bundle
elif [ $PLATFORM == "ILLUMINA" ];then
    export ML_MODEL=/resource/DNAscopeIlluminaWGS2.0.bundle
else
    echo $PLATFORM unknown
    exit 1
fi

if [ $KEEP_CLEAN == "not" ];then
    saferm $clean1 $clean2
fi

bed_arg=""
type_arg="WGS"
if [ $BED != "none" ];then
    bed_arg=" -b $BED"
    type_arg="WES"
fi

pcr_arg=""
if [ $PCR_FREE == "true" ];then
    pcr_arg=" --pcr_free "
fi

suffix_arg=""
if [ $SUFFIX == "bam" ];then
    suffix_arg=" --bam_format "
fi

[ -e dnascope.ok ]||dnascope
