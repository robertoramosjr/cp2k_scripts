#!/bin/bash
# Harpia job array for the Ag2SeO3 benchmark. One line of the list = "<abs_dir> <input.inp>".
# Usage: sbatch --array=0-35%3 [-p par16] run_array_harpia.sh <list_file>
# One node only on harpia (48 CPUs): %3 concurrent tasks of 16 CPUs each.
# RAM is 900 MB/CPU on harpia -> 16 CPUs = 14400M (--mem and --mem-per-cpu are exclusive).
#SBATCH -J ag2seo3
#SBATCH -p par16
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task=2
#SBATCH --mem=14400M
#SBATCH -t 3-00:00:00
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
set -u
LIST=$1
read -r DIR INP < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")
[ -z "${DIR:-}" ] && { echo "no line $SLURM_ARRAY_TASK_ID in $LIST"; exit 1; }
module load modmat/cp2k/2026.2-intel
ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_STACKSIZE=512M
export OMP_PLACES=cores OMP_PROC_BIND=close MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export I_MPI_PIN=1 I_MPI_PIN_DOMAIN=omp I_MPI_PIN_ORDER=compact
export FI_PROVIDER=shm
export CP2K_DATA_DIR=/home/hpc/software/cp2k/2026.2-intel/share/cp2k/data
cd "$DIR" || exit 1
PROJ=${INP%.inp}
# restart-safe: skip finished runs; resume GEO_OPT from the .restart
if [ -f "$PROJ.out" ] && tail -n 400 "$PROJ.out" | grep -q -E "OPTIMIZATION COMPLETED|PROGRAM ENDED"; then
  grep -q "SCF run NOT converged" "$PROJ.out" || { echo "[skip] $DIR done"; exit 0; }
fi
[ -f "${PROJ}-1.restart" ] && INP="${PROJ}-1.restart"
echo "task=$SLURM_ARRAY_TASK_ID dir=$DIR inp=$INP ranks=$SLURM_NTASKS threads=$OMP_NUM_THREADS node=$(hostname)"
mpirun -n $SLURM_NTASKS cp2k.psmp -i "$INP" -o "$PROJ.out"
grep -E "Total number of message passing|OPTIMIZATION COMPLETED|SCF run NOT converged|PROGRAM ENDED" "$PROJ.out" | tail -5
