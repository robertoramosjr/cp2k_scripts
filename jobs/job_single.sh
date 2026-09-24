#!/bin/bash
#SBATCH -J cp2k-run
#SBATCH -p short
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2194M
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err
#
# One run directory (steps 03, 04, 05, or one slab from 07):
#   sbatch ~/work_cp2k/jobs/job_single.sh RUN_DIR
# Restart-safe: resubmit the same command after a walltime kill and it
# continues from <project>-1.restart.

source "$(dirname "$(readlink -f "$0")")/common.sh" 2>/dev/null || source ~/work_cp2k/jobs/common.sh

run_dir "$(readlink -f "${1:-$SLURM_SUBMIT_DIR}")"
