# jobs/common.sh -- sourced by every job script. Two clusters:
#   access  (GridUNESP): 56 cores / 126 GB per node, partitions short (1 d) / medium (7 d) / long (30 d);
#                        module cp2k/2026.1 (MPICH) launched with Hydra mpiexec. Set up below.
#   coaraci (Unicamp):   everything cluster-specific is in common_coaraci.sh (whole-node rule:
#                        tasks x threads must be 48; local CP2K 2026.2 ssmp). Submit with sbatch_coaraci.sh.
# Chosen by CP2K_CLUSTER=access|coaraci, else by partition name (see below).

ulimit -s unlimited

# Cluster choice. Preferred: explicit CP2K_CLUSTER (sbatch_coaraci.sh exports
# CP2K_CLUSTER=coaraci, like the coaraci TiO2 job.sh files hardcode their
# environment). Fallback for plain sbatch: the partition name (coaraci has
# par48-x/par480-x/teste/...; access has short/medium/long). Host names are NOT
# reliable: on coaraci compute nodes neither SLURM_CLUSTER_NAME nor the host
# name said "coaraci" inside a batch job (launch test 1854499 ran as access).
if [ -z "${CP2K_CLUSTER:-}" ]; then
    case "${SLURM_JOB_PARTITION:-}" in
        par48-*|par480-*|teste|paralela|fat|serial|clauber|gpu-*) CP2K_CLUSTER=coaraci ;;
        short|medium|long)                                          CP2K_CLUSTER=access ;;
        *) case "$(hostname -f 2>/dev/null) $(getent hosts "$(hostname)" 2>/dev/null)" in
               *coaraci*) CP2K_CLUSTER=coaraci ;;
               *)         CP2K_CLUSTER=access ;;
           esac ;;
    esac
fi
echo "[INFO] common.sh: cluster=$CP2K_CLUSTER (partition ${SLURM_JOB_PARTITION:-none}, host $(hostname))" >&2

if [ "$CP2K_CLUSTER" = coaraci ]; then
    source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/common_coaraci.sh" || exit 2
else
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    module purge
    module load cp2k/2026.1
    export CP2K_DATA_DIR="$CP2K_DATA"        # module sets CP2K_DATA; cp2k.psmp reads CP2K_DATA_DIR
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
fi

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
    echo "[INFO] $(date '+%F %T') $dir: $CP2K_CLUSTER $CP2K_EXE, ${SLURM_NTASKS:-1} ranks x $OMP_NUM_THREADS threads, input $input"
    run_cp2k "$input" "$out"
    if is_done "$inp" "$out"; then
        echo "OK" > status.txt
    elif [ -s "$rst" ]; then
        echo "IN_PROGRESS (resubmit to resume from $rst)" > status.txt; return 1
    else
        echo "FAILED (see $out)" > status.txt; return 1
    fi
}
