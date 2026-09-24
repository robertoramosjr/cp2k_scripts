#!/bin/bash
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -c 2
#SBATCH --partition=long
#SBATCH --mail-user=roberto.ramos@unesp.br
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o job_etot.out
#SBATCH -e job_etot.error
#SBATCH -J TiO2-etot
#
# Single production run (no array): 28 MPI ranks x 2 OpenMP threads = 56
# cores, filling one full node -- follows the 2-threads/rank starting point
# from cp2k.org's own FAQ recommendation. Since this is the only job running
# in this allocation (unlike the scan arrays), there's no throughput reason
# to hold back cores here -- give it the whole node.
#
# If your system ends up much bigger later (~80-120 atoms, production
# research systems), re-benchmark: parallel efficiency generally improves
# with system size, so you may want -N 2 (or more) at that point. For this
# small anatase test system, 1 node is already generous.

ulimit -s unlimited
ulimit -l unlimited
ulimit -a

module purge
module load miniconda
module load cp2k/2026.1

# ==========================================================
# ⚙️  RUNTIME CONFIG
# ==========================================================
CP2K_EXE="cp2k.psmp"
export CP2K_DATA_DIR="$CP2K_DATA"   

# CP2K 2026.1 is linked against its own toolchain MPICH 4.3.2 (not Intel MPI):
# plain `srun` (no --mpi plugin) started N independent 1-rank copies of the
# same run, all writing the same .out ("message passing processes 1"). MPICH's
# own mpiexec (Hydra, SLURM-aware) launches one real N-rank job.
check_mpi_ranks() {
    local n
    n=$(grep -m1 "Total number of message passing processes" "$1" 2>/dev/null | awk '{print $NF}')
    if [ "$n" != "$SLURM_NTASKS" ]; then
        echo "[WARNING] $1 reports ${n:-?} MPI ranks but SLURM_NTASKS=$SLURM_NTASKS -- launcher/MPI mismatch." >&2
    fi
}
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK #limits the number of threads used 
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores #avoids hyper-threading
export OMP_PROC_BIND=close #putz processes in the nearest cores available# module only sets CP2K_DATA; cp2k.psmp reads CP2K_DATA_DIR

# CP2K prints this banner when GEO_OPT/CELL_OPT converges normally.
# NOTE: verify against your cluster's actual output and adjust if needed.
CONVERGED_PATTERN="OPTIMIZATION COMPLETED"

# ==========================================================
# 📂 DIRECTORY MANAGEMENT
# ==========================================================
# Allows execution like: sbatch ~/scripts/job_etot.sh /path/to/etot_dir
# Expects exactly one .inp in that directory (from input_generator.py
# --run-type GEO_OPT).
RUN_DIR=${1:-$SLURM_SUBMIT_DIR}
cd "$RUN_DIR" || { echo "[ERROR] Directory not found: $RUN_DIR"; exit 1; }

orig_inp=$(ls -- *.inp 2>/dev/null | head -1)
if [ -z "$orig_inp" ]; then
    echo "[ERROR] No .inp file found in $RUN_DIR."
    exit 1
fi

project_name="${orig_inp%.inp}"
out_file="${orig_inp}.out"
restart_file="${project_name}-1.restart"

# --- Already converged from a previous submission? ---
if [ -s "$out_file" ] && grep -q "$CONVERGED_PATTERN" "$out_file"; then
    echo "[SUCCESS] $RUN_DIR already converged. Nothing to do."
    echo "OK (recovered)" > status.txt
    exit 0
fi

# --- Restart protocol (CP2K equivalent of cp CONTCAR POSCAR) ---
if [ -s "$restart_file" ]; then
    echo "[INFO] Found restart checkpoint. Resuming from $restart_file."
    run_input="$restart_file"
else
    echo "[INFO] No restart checkpoint found. Starting fresh."
    run_input="$orig_inp"
fi

echo "[INFO] Running Etot (GEO_OPT): $RUN_DIR (input: $run_input)"
{
    echo "-------------------------------------------"
    echo "  CP2K Etot (GEO_OPT) -- $RUN_DIR"
    echo "-------------------------------------------"
    echo " This job is running in dir: $(pwd)"
    echo " Using input: $run_input"
    echo " Ranks x threads: ${SLURM_NTASKS} x ${OMP_NUM_THREADS}"
    echo " Starting: $(date +%Hh-%Mmin-%Ss--%Y/%m/%d)"
    echo "-------------------------------------------"
} > step.log

mpiexec -n "$SLURM_NTASKS" "$CP2K_EXE" -i "$run_input" -o "$out_file"

check_mpi_ranks "$out_file"
{
    echo "-------------------------------------------"
    echo " End: $(date +%Hh-%Mmin-%Ss--%Y/%m/%d)"
    echo "-------------------------------------------"
} >> step.log

if grep -q "$CONVERGED_PATTERN" "$out_file"; then
    echo "[INFO] Converged."
    echo "OK" > status.txt
elif [ -s "$restart_file" ]; then
    echo "[INFO] Did not converge yet but wrote a restart checkpoint -- "
    echo "       resubmit this SAME job to continue from $restart_file."
    echo "IN_PROGRESS (restart checkpoint available)" > status.txt
    exit 1
else
    echo "[WARNING] Did not converge and produced no restart checkpoint -- check $out_file."
    echo "FAILED" > status.txt
    exit 1
fi

echo "[INFO] Etot finished. Relaxed structure is in $restart_file (or extract via extract_relaxed_structure.py)."