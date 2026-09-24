# jobs/common.sh -- sourced by every job script. Cluster: GridUNESP (access),
# 56 cores / 126 GB per node, partitions short (1 d) / medium (7 d) / long (30 d).

ulimit -s unlimited

module purge
module load cp2k/2026.1
export CP2K_DATA_DIR="$CP2K_DATA"        # module sets CP2K_DATA; cp2k.psmp reads CP2K_DATA_DIR
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
CP2K_EXE=cp2k.psmp

# cp2k/2026.1 is linked against its own toolchain MPICH 4.3.2. Plain `srun`
# (no PMI plugin) starts N independent 1-rank copies of the same run, all
# writing one .out ("message passing processes 1") -- that is how every
# legacy TiO2_smoke run went. MPICH's Hydra mpiexec is SLURM-aware and
# starts one N-rank job; -bind-to core:T pins each rank to its T OpenMP cores.
run_cp2k() {  # run_cp2k <input> <output>
    mpiexec -n "$SLURM_NTASKS" -bind-to "core:${OMP_NUM_THREADS}" "$CP2K_EXE" -i "$1" -o "$2"
    local n
    n=$(grep -m1 "Total number of message passing processes" "$2" 2>/dev/null | awk '{print $NF}')
    if [ "$n" != "$SLURM_NTASKS" ]; then
        echo "[WARNING] $2 reports ${n:-?} MPI ranks, SLURM_NTASKS=$SLURM_NTASKS: launcher/MPI mismatch." >&2
    fi
}

run_type_of() { grep -m1 -E "^\s*RUN_TYPE" "$1" | awk '{print $2}'; }

# Success means converged, never just "PROGRAM ENDED AT": CP2K 2026.1 aborts
# on an unconverged SCF unless IGNORE_CONVERGENCE_FAILURE is set (grid scan),
# and an optimizer that hits MAX_ITER still ends normally.
is_done() {  # is_done <input> <output>
    local rt; rt=$(run_type_of "$1")
    [ -s "$2" ] || return 1
    case "$rt" in
        GEO_OPT|CELL_OPT) grep -q "OPTIMIZATION COMPLETED" "$2" ;;
        *) grep -q "PROGRAM ENDED AT" "$2" && grep -q "ENERGY| Total FORCE_EVAL" "$2" ;;
    esac
}

# Runs one directory: resumes GEO/CELL_OPT from <project>-1.restart (the
# CP2K "cp CONTCAR POSCAR"), writes status.txt.
run_dir() {
    local dir=$1 inp out project rst input
    cd "$dir" || { echo "[ERROR] no dir $dir"; return 1; }
    inp=$(ls -- *.inp 2>/dev/null | head -1)
    [ -n "$inp" ] || { echo "NO_INPUT" > status.txt; return 1; }
    project=${inp%.inp}; out="${inp}.out"; rst="${project}-1.restart"
    if is_done "$inp" "$out"; then echo "OK (already done)" > status.txt; return 0; fi
    input=$inp
    if [ -s "$rst" ]; then
        [ -s "$out" ] && mv "$out" "${out}.$(date +%s)"
        input=$rst
        echo "[INFO] resuming from $rst"
    fi
    echo "[INFO] $(date '+%F %T') $dir: $SLURM_NTASKS ranks x $OMP_NUM_THREADS threads, input $input"
    run_cp2k "$input" "$out"
    if is_done "$inp" "$out"; then
        echo "OK" > status.txt
    elif [ -s "$rst" ]; then
        echo "IN_PROGRESS (resubmit to resume from $rst)" > status.txt; return 1
    else
        echo "FAILED (see $out)" > status.txt; return 1
    fi
}
