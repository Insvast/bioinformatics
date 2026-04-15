#!/bin/bash
#SBATCH --job-name=WL-WES-trio-test-HG002
#SBATCH --partition=deimos
#SBATCH --cpus-per-task=64
#SBATCH -o %x-%j.log
#SBATCH -e %x-%j.log

set -xeuo pipefail

soft=/miniconda3/envs/myenv
export PATH=/apps/pkg/bin:$soft/bin:$PATH
export LD_LIBRARY_PATH=/apps/pkg/bin:$soft/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=/apps/pkg/bin:$soft/lib:${LIBRARY_PATH:-}
SENTIEON_INSTALL_DIR=/APP/u22/x86_com/sentieon/202503.03/sentieon-genomics-202503.03
export PATH=$SENTIEON_INSTALL_DIR/bin:$PATH
export PATH=/apps/pkg/bin/run-t1k:$PATH
export PATH=/apps/pkg/bin/samtools:$PATH
export PATH=/miniconda3/bin/sentieon-cli:$PATH
export LD_PRELOAD=/miniconda3/envs/anasys/lib/libjemalloc.so.2



SAMPLEID=HG002
WORKDIR=/project/WES-trio-test/result/$SAMPLEID
FASTQ_1=/project/WES-trio-test/data/HG002/HG002.WES_R1.fq.gz
FASTQ_2=/project/WES-trio-test/data/HG002/HG002.WES_R2.fq.gz
REFERENCE=/project/WES-trio-test/refseq/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
BED=/project/WES-trio-test/refseq/hg38.refGene.exon.bed
PCR_FREE=" --pcr-free "
# PLATFORM="ILLUMINA"  # "DNBSEQ"
PLATFORM="DNBSEQ"    #
export ML_MODEL=/CHSNP/app/resource/DNAscopeMGIWES2.1.bundle
# export ML_MODEL=/CHSNP/app/resource/DNAscopeIlluminaWES2.0.bundle
mkdir -p $WORKDIR
cd $WORKDIR

# -d dbsnp
sentieon-cli dnascope -r $REFERENCE --r1-fastq $FASTQ_1 --r2-fastq $FASTQ_2 \
  --readgroups "@RG\tID:$SAMPLEID\tSM:$SAMPLEID\tPL:$PLATFORM" \
  -m $ML_MODEL --duplicate-marking markdup \
  --assay WES --gvcf --skip-svs --skip-multiqc $PCR_FREE --bam_format $SAMPLEID.vcf.gz