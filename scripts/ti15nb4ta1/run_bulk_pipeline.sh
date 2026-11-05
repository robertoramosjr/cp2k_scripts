#!/bin/bash
# Harpia par48b, 48 CPUs (QoS harpia_48b: cpu=48, mem=48G/job; DefMemPerCPU 900 -> 43200M).
# 24 MPI ranks x 2 threads per full allocation (1 rank ~1 GiB; 48 ranks would exceed 43.2 GB).
# Runs bulk_pipeline.py (cutoff -> relcutoff -> kmesh -> cellopt); restart-safe: resubmit the same command.
# Usage: sbatch run_bulk_pipeline.sh [extra bulk_pipeline.py args]
#SBATCH -J tnt_bulk
#SBATCH -p par48b
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --cpus-per-task=1
#SBATCH --mem=43200M
#SBATCH -t 3-00:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -u
module load modmat/cp2k/2026.2-intel
ulimit -s unlimited
export CP2K_DATA_DIR=/home/hpc/software/cp2k/2026.2-intel/share/cp2k/data
export OMP_NUM_THREADS=2
cd /home/rramos/work_cp2k/Ti15Nb4Ta1 || exit 1
echo "job=$SLURM_JOB_ID cpus=$(nproc) node=$(hostname) $(date)"
python3 -u bulk_pipeline.py --stage all --ncpus "$SLURM_NTASKS" "$@"
