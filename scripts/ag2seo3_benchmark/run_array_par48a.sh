#!/bin/bash
# Harpia par48a array for the Ag2SeO3 k-mesh CELL_OPT benchmark (4 functionals).
# Half of the QoS (48 CPUs) per task = 24 CPUs = 12 MPI ranks x 2 threads; run 2 at a time (%2).
# RAM 900 MB/CPU on harpia -> 24 CPUs = 21600M (--mem and --mem-per-cpu are exclusive).
# One list line = "<abs_dir> <input.inp>".  Usage: sbatch --array=0-N%2 run_array_par48a.sh <list_file>
#SBATCH -J ag_kconv
#SBATCH -p par48a
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --cpus-per-task=2
#SBATCH --mem=21600M
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
# restart-safe: skip finished runs; resume CELL_OPT from the .restart
if [ -f "$PROJ.out" ] && tail -n 400 "$PROJ.out" | grep -q -E "OPTIMIZATION COMPLETED|PROGRAM ENDED"; then
  grep -q "SCF run NOT converged" "$PROJ.out" || { echo "[skip] $DIR done"; exit 0; }
fi
[ -f "${PROJ}-1.restart" ] && INP="${PROJ}-1.restart"
echo "task=$SLURM_ARRAY_TASK_ID dir=$DIR inp=$INP ranks=$SLURM_NTASKS threads=$OMP_NUM_THREADS node=$(hostname)"
mpirun -n $SLURM_NTASKS cp2k.psmp -i "$INP" -o "$PROJ.out"
grep -E "Total number of message passing|OPTIMIZATION COMPLETED|SCF run NOT converged|PROGRAM ENDED" "$PROJ.out" | tail -5
