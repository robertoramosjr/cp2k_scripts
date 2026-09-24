#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --partition=short
#SBATCH --mail-user=roberto.ramos@unesp.br
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o job_kmeshstress_%A_%a.out
#SBATCH -e job_kmeshstress_%A_%a.error
#SBATCH -J TiO2-kmesh-stress
#
# NOTE: --array is intentionally NOT set here (see job_cutoff_convergence_array.sh
# for the same rationale). Pass it at submission time:
#
#   N=$(find kmesh_stress_convergence -mindepth 1 -maxdepth 1 -type d | wc -l)
#   sbatch --array=0-$((N-1)) ~/scripts/job_kmesh_stress_convergence_array.sh \
#       $(pwd)/kmesh_stress_convergence
#
# Resource sizing (8 tasks x 2 cpus = 16 cores/array-task) targets the CURRENT
# small anatase test system. For production systems (~80-120 atoms), do a
# core-scaling test and bump --ntasks-per-node/--cpus-per-task toward a full
# node (or more) per array task -- see job_cutoff_convergence_array.sh header.
#
# Each CELL_OPT here is more expensive than a cutoff single-point (full
# relaxation), so this job requests a longer default walltime (24h) than the
# cutoff array job. If a task times out mid-relaxation, just resubmit the
# SAME array index -- the restart protocol below picks up from the .restart
# checkpoint automatically (see MOTION/PRINT/RESTART in input_generator.py).

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

# CP2K prints this banner when GEO_OPT/CELL_OPT converges normally -- but an
# ENERGY_FORCE (single-point) input never runs an optimizer at all, so this
# banner would never appear even on a perfectly good run. The right success
# criterion depends on RUN_TYPE, which is detected per-input below.
OPT_CONVERGED_PATTERN="OPTIMIZATION COMPLETED"
SCF_CONVERGED_PATTERN="outer SCF loop converged"
ENDED_PATTERN="PROGRAM ENDED AT"

# ==========================================================
# 📂 DIRECTORY MANAGEMENT (For Central Script Library Use)
# ==========================================================
SCAN_DIR=${1:-$SLURM_SUBMIT_DIR/kmesh_stress_convergence}

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
# 🚀 SINGLE K-MESH POINT (this array task's assigned directory)
# ==========================================================
orig_inp=$(ls -- *.inp 2>/dev/null | head -1)
if [ -z "$orig_inp" ]; then
    echo "[WARNING] No .inp file found in $run_dir."
    echo "NO_INPUT" > status.txt
    exit 1
fi

project_name="${orig_inp%.inp}"
out_file="${orig_inp}.out"
restart_file="${project_name}-1.restart"

# RUN_TYPE decides what "converged" even means: GEO_OPT/CELL_OPT run an
# optimizer and print an OPTIMIZATION COMPLETED banner when it converges,
# but ENERGY_FORCE/ENERGY are single-point -- no optimizer ever runs, so
# that banner never appears and SCF convergence is the only thing to check.
run_type=$(awk '/^[[:space:]]*RUN_TYPE[[:space:]]/{print $2; exit}' "$orig_inp")
is_converged() {
    local f="$1"
    if [ "$run_type" = "CELL_OPT" ] || [ "$run_type" = "GEO_OPT" ]; then
        grep -q "$OPT_CONVERGED_PATTERN" "$f"
    else
        grep -q "$ENDED_PATTERN" "$f" && grep -q "$SCF_CONVERGED_PATTERN" "$f"
    fi
}

# --- Already converged from a previous submission of this same array index? ---
if [ -s "$out_file" ] && is_converged "$out_file"; then
    echo "[SUCCESS] $run_dir already converged. Skipping."
    echo "OK (recovered)" > status.txt
    exit 0
fi

# --- Restart protocol (CP2K equivalent of cp CONTCAR POSCAR): -----------
# If a .restart file exists from an interrupted run, it IS a full valid
# CP2K input (last converged cell + coordinates + step counter) -- just
# run it directly instead of the original .inp.
if [ -s "$restart_file" ]; then
    echo "[INFO] Found restart checkpoint for $run_dir. Resuming from $restart_file."
    run_input="$restart_file"
else
    echo "[INFO] No restart checkpoint found for $run_dir. Starting fresh."
    run_input="$orig_inp"
fi

echo "[INFO] Running k-mesh stress point: $run_dir (array task $SLURM_ARRAY_TASK_ID, input: $run_input)"
{
    echo "-------------------------------------------"
    echo "  CP2K K-mesh Stress Convergence -- $run_dir (array task $SLURM_ARRAY_TASK_ID)"
    echo "-------------------------------------------"
    echo " This job is running in dir: $(pwd)"
    echo " Using input: $run_input"
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

if is_converged "$out_file"; then
    echo "[INFO] $run_dir converged."
    echo "OK" > status.txt
elif [ -s "$restart_file" ]; then
    echo "[INFO] $run_dir did not converge yet but wrote a restart checkpoint -- "
    echo "       resubmit this SAME array index to continue from $restart_file."
    echo "IN_PROGRESS (restart checkpoint available)" > status.txt
    exit 1
else
    echo "[WARNING] $run_dir did not converge and produced no restart checkpoint -- check $out_file."
    echo "FAILED" > status.txt
    exit 1
fi

echo "[INFO] Array task $SLURM_ARRAY_TASK_ID ($run_dir) done."
echo "[INFO] Once ALL array tasks show status.txt=OK, parse with:"
echo "  python ~/scripts/parse_kmesh_stress_convergence.py --scan-dir $SCAN_DIR --n-atoms <N_ATOMS>"