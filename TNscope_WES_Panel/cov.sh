#!/bin/bash
#SBATCH --job-name=WL-TumorNormal-cov
#SBATCH --partition=deimos
#SBATCH -n 64
#SBATCH -o %x-%j.log
#SBATCH -e %x-%j.log

echo ">>> 脚本开始时间：$(date "+%Y-%m-%d %H:%M:%S")"

bed=AIExome_Human_Exome_Panel_V3_Inherit-T600V1G-hg38.target.bed
bam=tumor_deduped.bam
out=tumor
mkdir -p $out
bamdst -p $bed -o $out $bam

bed=AIExome_Human_Exome_Panel_V3_Inherit-T600V1G-hg38.target.bed
bam=normal_deduped.bam
out=normal
mkdir -p $out
bamdst -p $bed -o $out $bam

echo ">>> 脚本结束时间：$(date "+%Y-%m-%d %H:%M:%S")"
