#!/bin/bash
# Production Etot single point (RUN_TYPE ENERGY, 6x5x5, relaxed 6x5x5 CELL_OPT geometry, no CELL_REF).
# Harpia par48b, 48 CPUs, 24 MPI ranks x 2 threads (same launch as bulk_pipeline.Pool; 1 GiB/rank measured).
# Success = PROGRAM ENDED + energy in the LAST output segment and no "SCF run NOT converged".
# Usage: sbatch run_etot.sh [rundir]   (default bulk/05_etot_k6-5-5; input written by make_etot.py)
#SBATCH -J tnt_etot
#SBATCH -p par48b
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --cpus-per-task=1
#SBATCH --mem=43200M
#SBATCH -t 0-12:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -u
module load modmat/cp2k/2026.2-intel
ulimit -s unlimited
export CP2K_DATA_DIR=/home/hpc/software/cp2k/2026.2-intel/share/cp2k/data
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OMP_PROC_BIND=false I_MPI_PIN=0 I_MPI_HYDRA_BOOTSTRAP=fork FI_PROVIDER=shm OMP_STACKSIZE=512M
RUN=/home/rramos/work_cp2k/Ti15Nb4Ta1/${1:-bulk/05_etot_k6-5-5}
P=TiNbTa_etot
cd "$RUN" || exit 1
echo "job=$SLURM_JOB_ID node=$(hostname) $(date) run=$RUN"
mpirun -n $((SLURM_NTASKS / OMP_NUM_THREADS)) cp2k.psmp -i $P.inp -o $P.out > $P.stdout 2>&1
seg=$(awk '/PROGRAM STARTED AT/{s=NR} END{print s}' $P.out)
tail -n +$seg $P.out > last_segment.out
grep -q "PROGRAM ENDED AT" last_segment.out || { echo "FAIL: no PROGRAM ENDED"; exit 1; }
grep -q "SCF run NOT converged" last_segment.out && { echo "FAIL: SCF NOT converged"; exit 1; }
grep -q "Total number of message passing processes *$((SLURM_NTASKS / OMP_NUM_THREADS))$" last_segment.out || echo "WARN: rank count mismatch"
E=$(grep "ENERGY| Total FORCE_EVAL" last_segment.out | tail -1 | awk '{print $NF}')
[ -n "$E" ] || { echo "FAIL: no energy"; exit 1; }
python3 - "$E" <<'PY'
import sys, json, re
E = float(sys.argv[1]); m = json.load(open("etot_meta.json")); t = open("last_segment.out").read()
pr = re.findall(r"STRESS\|\s+1/3 Trace\s+([-\d.Ee+]+)", t)
nscf = re.findall(r"SCF run converged in\s+(\d+) steps", t)
m.update(E_ha=E, E_ev=E * 27.211386245988, E_ev_per_atom=E * 27.211386245988 / m["natoms"],
         pressure_trace_bar=float(pr[-1]) if pr else None, scf_steps=int(nscf[-1]) if nscf else None)
json.dump(m, open("etot.json", "w"), indent=1); print(json.dumps(m, indent=1))
PY
echo ALL_OK
