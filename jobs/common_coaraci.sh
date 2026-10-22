# jobs/common_coaraci.sh -- coaraci (Unicamp) environment, sourced by common.sh
# when SLURM_CLUSTER_NAME is coaraci* (or CP2K_CLUSTER=coaraci).
#
# Node: 48 cores (2 x 24, no SMT), 120400 MB. Slurm ALWAYS hands the job the
# whole node: cores not requested sit idle and are charged anyway. So every job
# must use all 48 cores: tasks per node x cpus per task = 48. A job that does
# not is stopped here, before CP2K starts.
#
#   1 task  x 48 threads -> local CP2K 2026.2 cp2k.ssmp (OpenMP only; the default)
#   N tasks x 48/N       -> cp2k/2025.1 module cp2k.psmp (OpenMPI 5 + pmix), srun --mpi=pmix
#                           (different CP2K version: do not mix with 2026.2 runs of one project)
#
# Easiest: submit with jobs/sbatch_coaraci.sh (adds par48-x, 3 days, 1 x 48, memory), e.g.
#   ~/work_cp2k/cp2k_scripts/jobs/sbatch_coaraci.sh -J name ~/work_cp2k/cp2k_scripts/jobs/job_sequential.sh LIST

CORES_PER_NODE=${CORES_PER_NODE:-48}

_tpn=${SLURM_NTASKS_PER_NODE:-$(( ${SLURM_NTASKS:-1} / ${SLURM_NNODES:-1} ))}
_cpt=${SLURM_CPUS_PER_TASK:-1}
if [ -n "${SLURM_JOB_ID:-}" ] && [ $(( _tpn * _cpt )) -ne "$CORES_PER_NODE" ]; then
    echo "[ERROR] coaraci allocates the whole node: ${_tpn} tasks/node x ${_cpt} cpus/task = $(( _tpn * _cpt )) cores," >&2
    echo "        not ${CORES_PER_NODE}; the rest would sit idle. Resubmit with tasks x cpus = ${CORES_PER_NODE}" >&2
    echo "        (e.g. --ntasks-per-node=1 --cpus-per-task=48, or jobs/sbatch_coaraci.sh)." >&2
    exit 2
fi

export OMP_NUM_THREADS=$_cpt
export OMP_PLACES=cores
export OMP_PROC_BIND=close

if [ "${SLURM_NTASKS:-1}" -eq 1 ]; then
    # Static prebuilt 2026.2 (env.sh sets PATH, CP2K_DATA_DIR, OMP_STACKSIZE). One process,
    # so the srun singleton trap below does not apply.
    source "${CP2K_ENV:-$HOME/software/cp2k-2026.2/env.sh}"
    CP2K_EXE=cp2k.ssmp
    run_cp2k() {  # run_cp2k <input> <output>
        "$CP2K_EXE" -i "$1" -o "$2"
    }
else
    source /etc/profile.d/lmod.sh 2>/dev/null || source /etc/profile.d/*lmod* 2>/dev/null
    module load cp2k/2025.1-gcc-12.2.0-fod7ldk
    export CP2K_DATA_DIR=${CP2K_DATA_DIR:-/opt/spack/opt/spack/linux-oracle8-zen2/gcc-12.2.0/cp2k-2025.1/share/cp2k/data}
    export SRUN_CPUS_PER_TASK=$_cpt     # Slurm >= 22.05: else every rank is confined to 1 core
    CP2K_EXE=cp2k.psmp
    echo "[INFO] ${SLURM_NTASKS} tasks: using cp2k/2025.1 psmp (NOT the 2026.2 ssmp)" >&2
    # Without --mpi=pmix every rank runs as an independent singleton copy.
    run_cp2k() {  # run_cp2k <input> <output>
        srun --mpi=pmix "$CP2K_EXE" -i "$1" -o "$2"
        local n
        n=$(grep -m1 "Total number of message passing processes" "$2" 2>/dev/null | awk '{print $NF}')
        if [ "$n" != "$SLURM_NTASKS" ]; then
            echo "[WARNING] $2 reports ${n:-?} MPI ranks, SLURM_NTASKS=$SLURM_NTASKS: launcher/MPI mismatch." >&2
        fi
    }
fi
unset _tpn _cpt
