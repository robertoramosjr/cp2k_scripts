#!/bin/bash
# Harpia (modmat/cp2k/2026.2-intel, Intel MPI). Usage: sbatch [opts] run_half.sh <dir> <inp>
# e.g. sbatch -p teste -t 00:30:00 jobs_harpia/run_half.sh /abs/dir/half_sh_trans_timing trans_sh_half.inp
#SBATCH -J azo_half
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -u
DIR=$1; INP=$2
module load modmat/cp2k/2026.2-intel
ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_STACKSIZE=512M
export OMP_PLACES=cores OMP_PROC_BIND=close MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export I_MPI_PIN=1 I_MPI_PIN_DOMAIN=omp I_MPI_PIN_ORDER=compact I_MPI_DEBUG=4
export FI_PROVIDER=shm
export CP2K_DATA_DIR=/home/hpc/software/cp2k/2026.2-intel/share/cp2k/data
cd "$DIR" || exit 1
# restart-safe: se o .restart existir, continua dele
PROJ=${INP%.inp}; RST=$(ls ${PROJ%_half}_half-1.restart 2>/dev/null | head -1)
[ -n "$RST" ] && INP=$RST
echo "ranks=$SLURM_NTASKS threads=$OMP_NUM_THREADS node=$(hostname)"
mpirun -n $SLURM_NTASKS cp2k.psmp -i "$INP" -o "${PROJ}.out"
grep -E "Total number of message passing|OPTIMIZATION COMPLETED|SCF run NOT converged|PROGRAM ENDED" "${PROJ}.out" | tail -5
