$ cat HG002.pangenome.sh
#!/bin/bash
#SBATCH --job-name=WL-NIST-HG002-hybird-PB
#SBATCH --partition=deimos
#SBATCH --cpus-per-task=64
#SBATCH -o %x-%j.log
#SBATCH -e %x-%j.log

set -xeuo pipefail


export SENTIEON_INSTALL_DIR=/APP/u22/x86_com/sentieon/202503.03/sentieon-genomics-202503.03
export LD_LIBRARY_PATH=$SENTIEON_INSTALL_DIR/lib
export SENTIEON_LICENSE=X.X.X.X:8091 # 需替换
export LD_PRELOAD=/miniconda3/envs/anasys/lib/libjemalloc.so.2
export PATH=$SENTIEON_INSTALL_DIR/bin:$PATH

export THREADS=$(nproc)
export SENTIEON_TMPDIR=/project/NIST/result.pangenome
export MODEL_BUNDLE=/CHSNP/app/resource/HybridIlluminaPacBio1.1.bundle
REFERENCE=/project/WES-trio-test/refseq/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
SAMPLEID=HG002_35X_ill_15X_PB

echo ">>> dnascope-hybrid 开始时间：$(date "+%Y-%m-%d %H:%M:%S")"


srbam=/project/NIST/result.illumina/raw/HG002_deduped.cram
lrbam=/project/NIST/result.pb/HG002_PB_10X/HG002_PB_10X_mm2_sorted_fq_0.cram

sentieon-cli dnascope-hybrid \
    -r $REFERENCE \
    --sr_aln $srbam \
    --lr_aln $lrbam \
    --rgsm RGSM \
    -m $MODEL_BUNDLE \
    -t $THREADS \
    -g \
    --skip_multiqc \
    ${SAMPLEID}.vcf.gz

md5sum ${SAMPLEID}.vcf.gz > ${SAMPLEID}.vcf.gz.md5

echo ">>> 结束时间：$(date "+%Y-%m-%d %H:%M:%S")"