#!/bin/bash
# ACCESS (cp2k/2026.1, MPICH) array for the Ag2SeO3 MD stage. UNTESTED: written from the
# cp2k-pipeline skill facts, not yet run on access. One list line = "<abs_dir> <input.inp>".
# Usage: sbatch --array=0-N%K [-N nodes] [-p long] run_array_access.sh <list_file>
# 192 atoms / ADDED_MOS 864: override resources at submit time (e.g. -N 2 --ntasks-per-node=14).
#SBATCH -J ag_md
#SBATCH -p long
#SBATCH -N 1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2194M
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
set -u
LIST=$1
read -r DIR INP < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")
[ -z "${DIR:-}" ] && { echo "no line $SLURM_ARRAY_TASK_ID in $LIST"; exit 1; }
module load cp2k/2026.1
export CP2K_DATA_DIR="$CP2K_DATA"
ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_STACKSIZE=512M
cd "$DIR" || exit 1
PROJ=${INP%.inp}
if [ -f "$PROJ.out" ] && tail -n 400 "$PROJ.out" | grep -q "PROGRAM ENDED"; then
  grep -q "SCF run NOT converged" "$PROJ.out" || { echo "[skip] $DIR done"; exit 0; }
fi
[ -f "${PROJ}-1.restart" ] && INP="${PROJ}-1.restart"   # MD restart (written every 50 steps)
echo "task=$SLURM_ARRAY_TASK_ID dir=$DIR inp=$INP ranks=$SLURM_NTASKS threads=$OMP_NUM_THREADS"
mpiexec -n $SLURM_NTASKS -bind-to core:$OMP_NUM_THREADS cp2k.psmp -i "$INP" -o "$PROJ.out"
grep -E "Total number of message passing|PROGRAM ENDED|SCF run NOT converged" "$PROJ.out" | tail -5
