#!/bin/bash
# Slab GEO_OPT, sequential over a list "DIR PROJECT NRANKS NTHREADS" (one per line; make_slab_inputs.py), restart-safe.
# Per slab: stage 1 PROJECT_s1.inp (ENERGY, no dipole, loose SCF) -> PROJECT_s1-RESTART.kp, copied to
# PROJECT-RESTART.kp; stage 2 PROJECT.inp (GEO_OPT, fixed-position dipole, SCF_GUESS RESTART).
# Finished slabs (OPTIMIZATION COMPLETED in the last .out segment, SCF converged) are skipped,
# unfinished ones resume from PROJECT-1.restart. Harpia par48a/b: 48 CPUs (NRANKS x NTHREADS = 48).
# Usage: sbatch -p par48a -J tnt_slabA run_slabs.sh slab_opt/list_A.txt
#SBATCH -J tnt_slab
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
export OMP_PROC_BIND=false I_MPI_PIN=0 I_MPI_HYDRA_BOOTSTRAP=fork FI_PROVIDER=shm OMP_STACKSIZE=512M
ROOT=/home/rramos/work_cp2k/Ti15Nb4Ta1
LIST=$ROOT/${1:?list file}
echo "job=$SLURM_JOB_ID node=$(hostname) $(date) list=$LIST"
done_ok() {  # $1 = .out
  [ -f "$1" ] || return 1
  local seg; seg=$(awk '/PROGRAM STARTED AT/{s=NR} END{print s+0}' "$1")
  tail -n +"$seg" "$1" > "$1.last"
  grep -q "OPTIMIZATION COMPLETED" "$1.last" && grep -q "PROGRAM ENDED AT" "$1.last" && ! grep -q "SCF run NOT converged" "$1.last"
}
while read -r DIR P NR NT; do
  [ -n "${DIR:-}" ] || continue
  cd "$ROOT/$DIR" || { echo "[FAIL] no dir $DIR"; continue; }
  if done_ok $P.out; then echo "[skip] $DIR done"; continue; fi
  export OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
  if [ ! -f $P-RESTART.kp ]; then
    if ! { grep -q "SCF run converged" ${P}_s1.out 2>/dev/null && [ -f ${P}_s1-RESTART.kp ]; }; then
      echo "[s1] $DIR ${NR}x${NT} $(date)"
      mpirun -n $NR cp2k.psmp -i ${P}_s1.inp -o ${P}_s1.out >> ${P}_s1.stdout 2>&1 < /dev/null
      grep -q "SCF run converged" ${P}_s1.out && [ -f ${P}_s1-RESTART.kp ] || { echo "[FAIL] $DIR stage 1"; continue; }
    fi
    cp ${P}_s1-RESTART.kp $P-RESTART.kp
  fi
  IN=$P.inp; [ -f $P-1.restart ] && IN=$P-1.restart
  echo "[start] $DIR input=$IN ${NR}x${NT} $(date)"
  t0=$(date +%s)
  mpirun -n $NR cp2k.psmp -i $IN -o $P.out >> $P.stdout 2>&1 < /dev/null  # mpirun would eat the list on stdin
  rc=$?
  grep -q "Total number of message passing processes *$NR$" $P.out || echo "[WARN] rank count mismatch in $DIR"
  if done_ok $P.out; then st=ok; else st=FAIL; fi
  echo "[$st] $DIR rc=$rc $(( $(date +%s) - t0 ))s $(grep 'ENERGY| Total FORCE_EVAL' $P.out.last 2>/dev/null | tail -1 | awk '{print $NF}')"
done < "$LIST"
echo "END $(date)"
