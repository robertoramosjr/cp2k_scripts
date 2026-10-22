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

for _c in "$(dirname "$(readlink -f "$0")")" "${CP2K_SCRIPTS:-/nonexistent}/jobs" ~/work_cp2k/jobs ~/work_cp2k/cp2k_scripts/jobs; do
    [ -s "$_c/common.sh" ] && grep -q "^run_dir()" "$_c/common.sh" && { source "$_c/common.sh"; break; }
done
declare -F run_dir >/dev/null || { echo "[ERROR] jobs/common.sh not found (set CP2K_SCRIPTS=<repo>)"; exit 2; }

run_dir "$(readlink -f "${1:-$SLURM_SUBMIT_DIR}")"
