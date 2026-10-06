#!/bin/bash
#SBATCH -J cp2k-seq
#SBATCH -p short
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2194M
#SBATCH -o seq_%j.out
#SBATCH -e seq_%j.err
#
# Several run directories one after another in ONE allocation (for clusters
# where arrays are unreliable or per-user limits make them pointless).
#   sbatch ~/work_cp2k/jobs/job_sequential.sh LIST_FILE
# LIST_FILE: one run directory per line, in the order to run (# comments ok).
# A failed/unfinished directory does not stop the rest. Restart-safe: after a
# walltime kill resubmit the same command; finished directories are skipped
# (status OK) and the interrupted one resumes from <project>-1.restart.
# coaraci: sbatch -p par48-x -t 3-00:00:00 --ntasks-per-node=1 --cpus-per-task=48 --mem-per-cpu=2400M ...

for _c in "$(dirname "$(readlink -f "$0")")" "${CP2K_SCRIPTS:-/nonexistent}/jobs" ~/work_cp2k/jobs ~/work_cp2k/cp2k_scripts/jobs; do
    [ -s "$_c/common.sh" ] && grep -q "^run_dir()" "$_c/common.sh" && { source "$_c/common.sh"; break; }
done
declare -F run_dir >/dev/null || { echo "[ERROR] jobs/common.sh not found (set CP2K_SCRIPTS=<repo>)"; exit 2; }

LIST=$(readlink -f "${1:?usage: sbatch job_sequential.sh LIST_FILE}")
mapfile -t DIRS < <(grep -v -E '^\s*(#|$)' "$LIST")
nfail=0
for d in "${DIRS[@]}"; do
    t0=$(date +%s)
    ( run_dir "$d" )                       # subshell: run_dir cd's into the directory
    rc=$?
    [ $rc -ne 0 ] && nfail=$((nfail + 1))
    echo "[SEQ] $(date '+%F %T') rc=$rc $(( ($(date +%s) - t0) / 60 )) min  $d: $(cat "$d/status.txt" 2>/dev/null)"
done
echo "[SEQ] done: ${#DIRS[@]} dirs, $nfail not OK"
[ $nfail -eq 0 ]
