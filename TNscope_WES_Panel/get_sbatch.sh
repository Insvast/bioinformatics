#!/bin/bash
#********************************************************************
# FileName: get_sbatch.sh
# Description: Generate one Slurm job for a paired tumor-normal
#              TNscope hybridization panel analysis.
#********************************************************************
set -euo pipefail

mkdir -p jobs jobs_completed


# Project parameters. Edit these for each new project, or export them before running.
PN=${PN:-TumorNormal-test}
dir=${dir:-/project/${PN}}
data_dir=${data_dir:-${dir}/data/panel}
result_dir=${result_dir:-${dir}/panel-result2}

tumor_sm=${tumor_sm:-tumor_sample}
normal_sm=${normal_sm:-normal_sample}

tumor_r1=${tumor_r1:-${data_dir}/tumor_panel_R1.fastq.gz}
tumor_r2=${tumor_r2:-${data_dir}/tumor_panel_R2.fastq.gz}
normal_r1=${normal_r1:-${data_dir}/normal_panel_R1.fastq.gz}
normal_r2=${normal_r2:-${data_dir}/normal_panel_R2.fastq.gz}

fasta=${fasta:-/reference/human/pangemone/hg38_ucsc.fa}
dbsnp=${dbsnp:-/reference/human/pangemone/Homo_sapiens_assembly38.dbsnp138.vcf.gz}
interval_file=${interval_file:-${data_dir}/extended_100bp.bed}
tnscope_filter=${tnscope_filter:-${result_dir}/tnscope_filter.py}

partition=${partition:-deimos}
threads=${threads:-64}

sentieon_root=${sentieon_root:-/APP/u22/x86_com/sentieon/202503.03/sentieon-genomics-202503.03}
sentieon_license=${sentieon_license:-"需要指定"}

script_dir=$(cd "$(dirname "$0")" && pwd)
tnscope_panel=${tnscope_panel:-${script_dir}/TNscpoe-WES-Panel.sh}
submit_dir=$(pwd)

mkdir -p "$result_dir"

for f in "$tnscope_panel" "$tumor_r1" "$tumor_r2" "$normal_r1" "$normal_r2" "$fasta" "$dbsnp" "$tnscope_filter"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: required file not found or empty: $f" >&2
    exit 1
  fi
done

if [ -n "$interval_file" ] && [ "$interval_file" != "NA" ] && [ ! -s "$interval_file" ]; then
  echo "ERROR: interval BED not found or empty: $interval_file" >&2
  exit 1
fi

job_name="${tumor_sm}-${normal_sm}.job"
job="jobs/${job_name}"

cat > "$job" <<EOF
#!/bin/bash
#SBATCH --job-name=WL-${PN}-TNscope-panel
#SBATCH --partition=${partition}
#SBATCH --cpus-per-task=${threads}
#SBATCH -o %x-%j.log
#SBATCH -e %x-%j.log

set -xeuo pipefail
cd "${submit_dir}"
echo "script start time: \$(date "+%Y-%m-%d %H:%M:%S")"

export SENTIEON_INSTALL_DIR=${sentieon_root}
export SENTIEON_LICENSE=${sentieon_license}
"${tnscope_panel}" "${PN}" "${result_dir}" "${tumor_sm}" "${tumor_r1}" "${tumor_r2}" "${normal_sm}" "${normal_r1}" "${normal_r2}" "${fasta}" "${dbsnp}" "${interval_file}" "${tnscope_filter}"

mkdir -p "${submit_dir}/jobs_completed"
mv "${submit_dir}/jobs/${job_name}" "${submit_dir}/jobs_completed/" || true
echo "script end time: \$(date "+%Y-%m-%d %H:%M:%S")"
EOF
chmod +x "$job"
echo "Generated $job"

cat > submit_tnscope_job.sh <<EOF
#!/bin/bash
set -euo pipefail

sbatch "jobs/${job_name}" | tee submitted.tnscope.jobs.txt
EOF
chmod +x submit_tnscope_job.sh
echo "Generated submit_tnscope_job.sh"
