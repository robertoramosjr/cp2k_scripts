#!/bin/bash
#SBATCH -t 5-00:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2194M
#SBATCH --partition=long
#SBATCH --mail-user=roberto.ramos@unesp.br
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o job_dos_bands.out
#SBATCH -e job_dos_bands.error
#SBATCH -J TiO2-dos-bands
#
# Two CP2K steps, ONE queue allocation:
#   1. verify-relax: quick GEO_OPT confirming the pymatgen-standardized cell
#      (used for the band path) is still at the energy minimum.
#   2. dos_bands: PDOS + band structure on the confirmed structure, with a
#      denser k-mesh than Etot used.
# Steps in between (extracting the relaxed CIF, building the DOS/bands
# input) are cheap Python and just re-run every time this job starts --
# only the two CP2K invocations have restart-skip logic.
#
# If walltime is hit before finishing everything, resubmit this SAME job:
# it detects how far it got and resumes from there (still one script, one
# config, no manual bookkeeping between steps -- but note this means "one
# submission" in the sense of "one job you manage", not an ironclad
# guarantee of finishing in a single queue slot if the run is unusually slow).

ulimit -s unlimited
ulimit -l unlimited
ulimit -a

module purge
module load miniconda
conda activate cp2k_env   # <-- ADJUST to whatever conda env has pymatgen/mp-api installed
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


CONVERGED_PATTERN="OPTIMIZATION COMPLETED"
ENERGY_DONE_PATTERN="PROGRAM ENDED AT"
# "PROGRAM ENDED AT" only means CP2K exited without crashing -- it does NOT
# abort on SCF non-convergence, it just prints whatever it has and moves on.
# Require the explicit outer-SCF convergence line too (STEP 4 is RUN_TYPE
# ENERGY, single-point).
SCF_CONVERGED_PATTERN="outer SCF loop converged"

# ==========================================================
# 🧾 USER CONFIG -- edit these per run (hardcoded on purpose: this is a
# pipeline/trigger job, not a reusable library script)
# ==========================================================
RELAXED_CIF="TiO2_relaxed.cif"        # structure coming out of Etot
PROJECT_NAME="TiO2_smoke"
CUTOFF=350
KPOINTS="9 4 4"                             # Etot's production mesh
VERIFY_FUNCTIONAL="PBE"    # STEP 1 (verify-relax): should match whatever level
                            # of theory produced RELAXED_CIF (Etot), so the
                            # confirmation is checking against the SAME PES it
                            # was relaxed on -- NOT necessarily the same as
                            # FUNCTIONAL below.
FUNCTIONAL="HSE06"         # STEP 3/4 (the actual DOS/bands SCF)
DOS_KPOINTS_MULTIPLIER=2.0
WFN_RESTART=""      # from Etot; leave empty "" to skip SCF_GUESS RESTART

# ==========================================================
# 📂 DIRECTORY MANAGEMENT
# ==========================================================
WORK_DIR=${1:-$SLURM_SUBMIT_DIR}
cd "$WORK_DIR" || { echo "[ERROR] Directory not found: $WORK_DIR"; exit 1; }

# ==========================================================
# 🚀 STEP 1: verify-relax on the pymatgen-standardized cell
# ==========================================================
echo "[STEP 1] Generating verify-relax input (functional: $VERIFY_FUNCTIONAL)..."
python ~/scripts/dos_bands_generator.py --cif "$RELAXED_CIF" --project-name "$PROJECT_NAME" \
    --cutoff "$CUTOFF" --kpoints $KPOINTS --functional "$VERIFY_FUNCTIONAL" --verify-relax

VERIFY_PROJECT="${PROJECT_NAME}_verify_relax"
VERIFY_INP="${VERIFY_PROJECT}.inp"
VERIFY_OUT="${VERIFY_INP}.out"
VERIFY_RESTART="${VERIFY_PROJECT}-1.restart"

if [ ! -s "$VERIFY_INP" ]; then
    echo "[ERROR] verify-relax input was not generated. Check the python call above."
    exit 1
fi

if [ -s "$VERIFY_OUT" ] && grep -q "$CONVERGED_PATTERN" "$VERIFY_OUT"; then
    echo "[STEP 1] Already converged. Skipping CP2K run."
else
    if [ -s "$VERIFY_RESTART" ]; then
        echo "[STEP 1] Resuming from restart checkpoint."
        verify_run_input="$VERIFY_RESTART"
    else
        echo "[STEP 1] Starting fresh."
        verify_run_input="$VERIFY_INP"
    fi

    echo "[STEP 1] Running verify-relax CP2K..."
    mpiexec -n "$SLURM_NTASKS" "$CP2K_EXE" -i "$verify_run_input" -o "$VERIFY_OUT"
    check_mpi_ranks "$VERIFY_OUT"
    if grep -q "$CONVERGED_PATTERN" "$VERIFY_OUT"; then
        echo "[STEP 1] Converged."
    elif [ -s "$VERIFY_RESTART" ]; then
        echo "[STEP 1] Did not converge yet, but wrote a restart checkpoint."
        echo "         Resubmit this SAME job to continue -- it will pick up from here."
        echo "STEP1_IN_PROGRESS" > status.txt
        exit 1
    else
        echo "[ERROR] verify-relax failed and produced no restart checkpoint -- check $VERIFY_OUT."
        echo "STEP1_FAILED" > status.txt
        exit 1
    fi
fi

# ==========================================================
# 📐 STEP 2: extract the confirmed relaxed structure
# ==========================================================
echo "[STEP 2] Extracting confirmed relaxed structure..."
CONFIRMED_CIF="${PROJECT_NAME}_standardized_relaxed.cif"
python ~/scripts/extract_relaxed_structure.py --restart "$VERIFY_RESTART" --output "$CONFIRMED_CIF"

if [ ! -s "$CONFIRMED_CIF" ]; then
    echo "[ERROR] Could not extract confirmed structure from $VERIFY_RESTART."
    echo "STEP2_FAILED" > status.txt
    exit 1
fi

# ==========================================================
# 🧮 STEP 3: build the DOS/bands input
# ==========================================================
echo "[STEP 3] Generating DOS/bands input (functional: $FUNCTIONAL)..."
WFN_FLAG=""
if [ -n "$WFN_RESTART" ] && [ -s "$WFN_RESTART" ]; then
    WFN_FLAG="--wfn-restart $WFN_RESTART"
fi

DOSBANDS_PROJECT="${PROJECT_NAME}_dos_bands"
python ~/scripts/dos_bands_generator.py --cif "$CONFIRMED_CIF" --project-name "$DOSBANDS_PROJECT" \
    --cutoff "$CUTOFF" --kpoints $KPOINTS --functional "$FUNCTIONAL" \
    --dos-kpoints-multiplier "$DOS_KPOINTS_MULTIPLIER" $WFN_FLAG

DOSBANDS_INP="${DOSBANDS_PROJECT}.inp"
DOSBANDS_OUT="${DOSBANDS_INP}.out"

if [ ! -s "$DOSBANDS_INP" ]; then
    echo "[ERROR] DOS/bands input was not generated. Check the python call above."
    echo "STEP3_FAILED" > status.txt
    exit 1
fi

# ==========================================================
# 🚀 STEP 4: run the DOS/bands calculation
# ==========================================================
if [ -s "$DOSBANDS_OUT" ] && grep -q "$ENERGY_DONE_PATTERN" "$DOSBANDS_OUT" \
   && grep -q "$SCF_CONVERGED_PATTERN" "$DOSBANDS_OUT"; then
    echo "[STEP 4] Already completed. Skipping CP2K run."
else
    echo "[STEP 4] Running DOS/bands CP2K..."
    mpiexec -n "$SLURM_NTASKS" "$CP2K_EXE" -i "$DOSBANDS_INP" -o "$DOSBANDS_OUT"
    check_mpi_ranks "$DOSBANDS_OUT"
fi

if ! grep -q "$ENERGY_DONE_PATTERN" "$DOSBANDS_OUT"; then
    echo "[WARNING] DOS/bands run did not report normal CP2K termination -- check $DOSBANDS_OUT."
    echo "STEP4_FAILED (crash/abort)" > status.txt
    exit 1
elif ! grep -q "$SCF_CONVERGED_PATTERN" "$DOSBANDS_OUT"; then
    echo "[WARNING] DOS/bands run ended but the SCF did NOT converge -- PDOS/bands are unreliable. Check $DOSBANDS_OUT."
    echo "STEP4_FAILED (SCF not converged)" > status.txt
    exit 1
else
    echo "[INFO] Pipeline finished. PDOS files: ${DOSBANDS_PROJECT}-k*-*.pdos, bands: ${DOSBANDS_PROJECT}.bs"
    echo "OK" > status.txt
fi