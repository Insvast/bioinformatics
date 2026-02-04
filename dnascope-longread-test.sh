#!/bin/bash
#SBATCH --job-name=WL-longread-test-XHTGS021
#SBATCH --partition=deimos
#SBATCH --cpus-per-task=64

# project/longread-test/code
# project/longread-test/refseq/GCF_000001405.40_GRCh38.p14_genomic.fna
# project/longread-test/data/XHTGS021/XHTGS021.no-signal.bam

source apps/myenv.sh
module load bcftools/1.21-gcc-11.4.0 samtools/1.20-gcc-11.4.0
export NT=64
export SENTIEON_INSTALL_DIR="/APP/u22/x86_com/sentieon/202503.01/sentieon-genomics-202503.01"
export MODEL_BUNDLE="/CHSNP/app/resource/DNAscopePacBio2.1.bundle"
export SENTIEON_EXEC="${SENTIEON_INSTALL_DIR}/bin/sentieon"

export SAMPLEID=XHTGS021
export PROJECTDIR="/project/longread-test"
# export REF="${PROJECTDIR}/refseq/GCF_000001405.40_GRCh38.p14_genomic.fna"
export REF="${PROJECTDIR}/refseq/hg38.fa"
export WORKDIR="${PROJECTDIR}/result/${SAMPLEID}"
mkdir -p $WORKDIR

cd $WORKDIR
# export BAM="${PROJECTDIR}/data/${SAMPLEID}/${SAMPLEID}.no-signal.bam"

export BAM_RAW="${PROJECTDIR}/data/${SAMPLEID}/${SAMPLEID}.no-signal.bam"

echo ">>> dnascope-longread 开始时间：$(date "+%Y-%m-%d %H:%M:%S")"
sentieon-cli dnascope-longread \
        -r ${REF} \
        -i ${BAM_RAW} \
        --align \
        -m ${MODEL_BUNDLE} \
        --tech HiFi \
        -t 60 \
        $WORKDIR/${SAMPLEID}.vcf.gz
md5sum $WORKDIR/${SAMPLEID}.vcf.gz > $WORKDIR/${SAMPLEID}.vcf.gz.md5


echo ">>> 结束时间：$(date "+%Y-%m-%d %H:%M:%S")"
