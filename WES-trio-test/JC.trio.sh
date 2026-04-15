#!/bin/bash
#SBATCH --job-name=WL-WES-trio-test
#SBATCH --partition=deimos
#SBATCH --cpus-per-task=64
#SBATCH -o %x-%j.log
#SBATCH -e %x-%j.log
set -xeuo pipefail

echo ">>> 脚本开始时间：$(date "+%Y-%m-%d %H:%M:%S")"

REFERENCE=/project/WES-trio-test/refseq/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
root=/APP/u22/x86_com/sentieon/202503.03/sentieon-genomics-202503.03
ML_MODEL=/CHSNP/app/resource/DNAscopeMGIWGS2.0.bundle


proband_name=HG002
father_name=HG003
mother_name=HG004
# 步骤2：家系联合基因型进行初轮联合变异检测

$root/bin/sentieon driver -r $REFERENCE --algo GVCFtyper \
  -v $proband_name.g.vcf.gz -v $father_name.g.vcf.gz -v $mother_name.g.vcf.gz \
  joint-call_pass1.vcf.gz

# 步骤3： de novo 突变检测

bcftools +trio-dnm2 -p $proband_name,$father_name,$mother_name -X GRCh38 \
  -o trio-dnm2.output.vcf.gz -Oz --use-NAIVE joint-call_pass1.vcf.gz

# 步骤4：变异位点分类

# 筛选孟德尔不兼容变异（可能是DNM）
bcftools view -i "FMT/DNM[0]==1" trio-dnm2.output.vcf.gz | \
  $root/bin/sentieon util vcfconvert - mendelian-incompatible.vcf.gz
# 筛选孟德尔兼容变异（符合遗传规律）
bcftools view -e "FMT/DNM[0]==1" trio-dnm2.output.vcf.gz | \
  $root/bin/sentieon util vcfconvert - mendelian-compatible.vcf.gz

# 步骤5：重召回孟德尔不兼容变异

$root/bin/sentieon driver -r $REFERENCE \
  -i ${proband_name}_deduped.bam -i ${father_name}_deduped.bam -i ${mother_name}_deduped.bam \
  --algo DNAscope --pcr_indel_model none --given mendelian-incompatible.vcf.gz \
  --model $ML_MODEL/dnascope.model mendelian-incompatible.recalled.vcf.gz

# 步骤6：筛选高置信度重召回结果

qual_thresh=20
bcftools view -e "QUAL<$qual_thresh" mendelian-incompatible.recalled.vcf.gz | \
  $root/bin/sentieon util vcfconvert - mendelian-incompatible.recalled.highconf.vcf.gz

# 步骤7：合并结果并进行第二次联合调用

bcftools concat --allow-overlaps mendelian-compatible.vcf.gz \
  mendelian-incompatible.recalled.highconf.vcf.gz | \
  bcftools sort - | $root/bin/sentieon util vcfconvert - joint-call_pass2.vcf.gz

# 步骤8：第二次检测 de novo 突变并建立索引

bcftools +trio-dnm2 -p $proband_name,$father_name,$mother_name -X GRCh38 \
  -o trio-dnm2.joint-call_pass2.vcf.gz -Oz --use-NAIVE joint-call_pass2.vcf.gz
$root/bin/sentieon util vcfindex trio-dnm2.joint-call_pass2.vcf.gz

echo ">>> 脚本结束时间：$(date "+%Y-%m-%d %H:%M:%S")"