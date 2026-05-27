#!/bin/bash
set -euo pipefail

sbatch "jobs/SRR7890852-SRR7890847.job" | tee submitted.tnscope.jobs.txt
