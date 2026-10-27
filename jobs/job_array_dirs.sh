#!/bin/bash
#SBATCH -J cp2k-arr
#SBATCH -p short
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2194M
#SBATCH -o job_%A_%a.out
#SBATCH -e job_%A_%a.err
#
# One run directory per array task, listed one per line in DIRLIST (task N = line N):
#   sbatch --array=1-2 ~/work_cp2k/jobs/job_array_dirs.sh DIRLIST
# Restart-safe like job_single.sh: resubmit the same command after a walltime kill.

for _c in "$(dirname "$(readlink -f "$0")")" "${CP2K_SCRIPTS:-/nonexistent}/jobs" ~/work_cp2k/jobs ~/work_cp2k/cp2k_scripts/jobs; do
    [ -s "$_c/common.sh" ] && grep -q "^run_dir()" "$_c/common.sh" && { source "$_c/common.sh"; break; }
done
declare -F run_dir >/dev/null || { echo "[ERROR] jobs/common.sh not found (set CP2K_SCRIPTS=<repo>)"; exit 2; }

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$(readlink -f "$1")")
[ -d "$DIR" ] || { echo "[ERROR] task $SLURM_ARRAY_TASK_ID: no directory '$DIR' in $1"; exit 2; }
run_dir "$(readlink -f "$DIR")"
