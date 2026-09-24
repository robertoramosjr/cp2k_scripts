#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --partition=short
#SBATCH --mail-user=roberto.ramos@unesp.br
#SBATCH --mail-type=END,FAIL
#SBATCH -o job_cutoffconv_%A_%a.out
#SBATCH -e job_cutoffconv_%A_%a.error
#SBATCH -J TiO2-cutoffconv
#
# NOTE: --array is intentionally NOT set here -- it depends on how many
# CUTOFF points exist in the scan, which is only known after
# cutoff_convergence.py has generated the directory tree. Pass it at
# submission time, e.g. for 9 scan points (indices 0-8):
#
#   N=$(find cutoff_convergence -mindepth 1 -maxdepth 1 -type d | wc -l)
#   sbatch --array=0-$((N-1)) ~/scripts/job_cutoff_convergence_array.sh \
#       $(pwd)/cutoff_convergence
#
# Resource sizing below (8 tasks x 2 cpus = 16 cores/array-task) targets the
# CURRENT small anatase test system (6-12 atoms), where per-calculation
# parallel efficiency saturates well before a full 56-core node -- so it's
# better to run several scan points concurrently than to over-allocate one.
# For production systems (~80-120 atoms), do a quick core-scaling test first
# and bump --ntasks-per-node / --cpus-per-task toward a full node (or more)
# per array task, since larger systems parallelize efficiently over more cores.

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

# "PROGRAM ENDED AT" only means CP2K exited without crashing -- it says
# nothing about whether the SCF actually converged (CP2K does NOT abort
# on SCF non-convergence, it just prints the unconverged energy and moves
# on). Require the explicit outer-SCF convergence line too, or a run that
# silently failed to converge would be reported as a normal success.
ENDED_PATTERN="PROGRAM ENDED AT"
SCF_CONVERGED_PATTERN="outer SCF loop converged"

# ==========================================================
# 📂 DIRECTORY MANAGEMENT (For Central Script Library Use)
# ==========================================================
SCAN_DIR=${1:-$SLURM_SUBMIT_DIR/cutoff_convergence}

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "[ERROR] This script must be submitted with --array=0-N (see header comment)."
    exit 1
fi

mapfile -t RUN_DIRS < <(find "$SCAN_DIR" -mindepth 1 -maxdepth 1 -type d | sort -V)

if [ "$SLURM_ARRAY_TASK_ID" -ge "${#RUN_DIRS[@]}" ]; then
    echo "[ERROR] Array index $SLURM_ARRAY_TASK_ID has no matching directory "
    echo "        (only ${#RUN_DIRS[@]} found in $SCAN_DIR). Check --array range."
    exit 1
fi

run_dir="${RUN_DIRS[$SLURM_ARRAY_TASK_ID]}"
cd "$run_dir" || { echo "[ERROR] Could not cd into $run_dir"; exit 1; }

# ==========================================================
# 🚀 SINGLE SCAN POINT (this array task's assigned directory)
# ==========================================================
inp_file=$(ls -- *.inp 2>/dev/null | head -1)
if [ -z "$inp_file" ]; then
    echo "[WARNING] No .inp file found in $run_dir."
    echo "NO_INPUT" > status.txt
    exit 1
fi

out_file="${inp_file}.out"

# --- Restart protocol: skip if already finished normally (e.g. array task
# was resubmitted after a partial previous run) ---
if [ -s "$out_file" ] && grep -q "$ENDED_PATTERN" "$out_file" && grep -q "$SCF_CONVERGED_PATTERN" "$out_file"; then
    echo "[SUCCESS] $run_dir already completed. Skipping."
    echo "OK (recovered)" > status.txt
    exit 0
fi

echo "[INFO] Running CUTOFF scan point: $run_dir (array task $SLURM_ARRAY_TASK_ID)"
{
    echo "-------------------------------------------"
    echo "  CP2K Cutoff Convergence -- $run_dir (array task $SLURM_ARRAY_TASK_ID)"
    echo "-------------------------------------------"
    echo " This job is running in dir: $(pwd)"
    echo " Starting: $(date +%Hh-%Mmin-%Ss--%Y/%m/%d)"
    echo "-------------------------------------------"
} > step.log

mpiexec -n "$SLURM_NTASKS" "$CP2K_EXE" -i "$inp_file" -o "$out_file"

check_mpi_ranks "$out_file"
{
    echo "-------------------------------------------"
    echo " End: $(date +%Hh-%Mmin-%Ss--%Y/%m/%d)"
    echo "-------------------------------------------"
} >> step.log

if ! grep -q "$ENDED_PATTERN" "$out_file"; then
    echo "[WARNING] $run_dir did not report normal CP2K termination -- check $out_file."
    echo "FAILED (crash/abort)" > status.txt
    exit 1
elif ! grep -q "$SCF_CONVERGED_PATTERN" "$out_file"; then
    echo "[WARNING] $run_dir ran to completion but the SCF did NOT converge -- check $out_file."
    echo "FAILED (SCF not converged)" > status.txt
    exit 1
else
    echo "[INFO] $run_dir finished normally and SCF converged."
    echo "OK" > status.txt
fi

echo "[INFO] Array task $SLURM_ARRAY_TASK_ID ($run_dir) done."
echo "[INFO] Once ALL array tasks finish, parse with:"
echo "  python ~/scripts/parse_cutoff_convergence.py --scan-dir $SCAN_DIR --n-atoms <N_ATOMS>"