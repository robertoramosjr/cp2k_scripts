#!/bin/bash
# probe_build.sh -- what does THIS CP2K build accept? (CP2K_PIPELINE.md, PART 3, items 1-2)
#
#   1. DOS/PDOS interface: legacy &DFT%PRINT%PDOS vs unified &DFT%PRINT%DOS(%PDOS)
#      -> decides 05_pdos.py --dos-interface
#   2. RI-HFX with k-points (&XC%HF%RI + &KPOINTS): syntax, and optionally a real run
#
# Safe mode (default, OK on the login node): schema dump (cp2k.psmp --xml) +
# `cp2k.psmp --check` on tiny inputs. Nothing is computed.
# --run: also RUNS the RI-HFX k-point test (Si, 2 atoms). RI-HFXk memory can
# reach tens of GiB per rank, so --run refuses to start outside a SLURM job:
#   sbatch -p short -N 1 --ntasks-per-node=4 --cpus-per-task=2 --mem-per-cpu=8G \
#          --wrap "bash ~/work_cp2k/scripts/probe_build.sh --run"
#
# Usage: bash probe_build.sh [--run] [--module cp2k/2026.1] [--workdir DIR]

set -u
MODULE=cp2k/2026.1
RUN=0
WORK=${TMPDIR:-/tmp}/cp2k_probe_$$
PY=${PYTHON:-$HOME/.conda/envs/cp2k_env/bin/python}
while [ $# -gt 0 ]; do
    case $1 in
        --run) RUN=1 ;;
        --module) MODULE=$2; shift ;;
        --workdir) WORK=$2; shift ;;
        *) echo "unknown option $1"; exit 2 ;;
    esac
    shift
done

module purge >/dev/null 2>&1
module load "$MODULE" >/dev/null 2>&1 || { echo "[ERROR] module $MODULE"; exit 1; }
export CP2K_DATA_DIR="$CP2K_DATA"
mkdir -p "$WORK" && cd "$WORK" || exit 1
echo "== probe of $MODULE in $WORK"

# ---------------------------------------------------------------- schema
cp2k.psmp --xml >/dev/null 2>&1
if [ ! -s cp2k_input.xml ]; then echo "[ERROR] cp2k.psmp --xml produced nothing"; exit 1; fi
"$PY" - <<'EOF'
import xml.etree.ElementTree as ET
root = ET.parse("cp2k_input.xml").getroot()
def sec(path):
    n = root
    for name in path.split("/"):
        nxt = [s for s in n.findall("SECTION") if s.findtext("NAME") == name]
        if not nxt:
            return None
        n = nxt[0]
    return n
def kws(n):
    return {k.findtext("NAME") for k in n.findall("KEYWORD")} if n is not None else set()
checks = {
    "legacy PDOS   FORCE_EVAL/DFT/PRINT/PDOS": "FORCE_EVAL/DFT/PRINT/PDOS",
    "unified DOS   FORCE_EVAL/DFT/PRINT/DOS": "FORCE_EVAL/DFT/PRINT/DOS",
    "unified PDOS  FORCE_EVAL/DFT/PRINT/DOS/PDOS": "FORCE_EVAL/DFT/PRINT/DOS/PDOS",
    "RI-HFX        FORCE_EVAL/DFT/XC/HF/RI": "FORCE_EVAL/DFT/XC/HF/RI",
}
for label, path in checks.items():
    print(f"  [schema] {label:<48} {'YES' if sec(path) is not None else 'no'}")
pd = sec("FORCE_EVAL/DFT/PRINT/PDOS")
if pd is not None:
    print(f"  [schema] legacy PDOS keywords: {sorted(kws(pd) & {'COMPONENTS', 'NLUMO', 'APPEND'})}")
ri = sec("FORCE_EVAL/DFT/XC/HF/RI")
print(f"  [schema] HF/RI keywords of interest: {sorted(kws(ri) & {'RI_FLAVOR', 'KP_RI_MEMORY_ESTIMATE', 'EPS_FILTER', 'RI_METRIC', 'MEMORY_CUT'})}")
EOF

# ---------------------------------------------------------------- tiny inputs
cell='    &CELL
      ABC 5.43 5.43 5.43
      ALPHA_BETA_GAMMA 90 90 90
    &END CELL
    &COORD
      SCALED
      Si 0.00 0.00 0.00
      Si 0.25 0.25 0.25
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-SR-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND'
head_dft='  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 200
    &END MGRID'
mk() {  # mk <name> <xc block> <extra DFT text> <kpoints yes/no>
    local kp=""
    [ "$4" = yes ] && kp='    &KPOINTS
      SCHEME MONKHORST-PACK 2 2 2
    &END KPOINTS'
    cat > "$1.inp" <<EOF2
&GLOBAL
  PROJECT $1
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
$head_dft
$2
    &SCF
      MAX_SCF 30
      ADDED_MOS 4
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
      &MIXING
        METHOD BROYDEN_MIXING
      &END MIXING
    &END SCF
$kp
$3
  &END DFT
  &SUBSYS
$cell
  &END SUBSYS
&END FORCE_EVAL
EOF2
}
PBE='    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC'
PBE0_RI='    &XC
      &XC_FUNCTIONAL
        &PBE
          SCALE_X 0.75
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &RI
        &END RI
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS 3.0
        &END INTERACTION_POTENTIAL
      &END HF
    &END XC'
mk pdos_legacy "$PBE" '    &PRINT
      &PDOS
        COMPONENTS
        NLUMO -1
      &END PDOS
    &END PRINT' yes
mk pdos_unified "$PBE" '    &PRINT
      &DOS
        &PDOS
          COMPONENTS
          NLUMO -1
        &END PDOS
      &END DOS
    &END PRINT' yes
mk rihfx_kp "$PBE0_RI" "" yes

for f in pdos_legacy pdos_unified rihfx_kp; do
    if cp2k.psmp --check "$f.inp" > "$f.check" 2>&1 && grep -q "SUCCESS" "$f.check"; then r=ACCEPTED; else r=REJECTED; fi
    echo "  [--check] $f.inp: $r$( [ $r = REJECTED ] && grep -m1 -E 'ABORT|unknown|not' $f.check | sed 's/^ */ -- /')"
done

# ---------------------------------------------------------------- optional run
if [ "$RUN" = 1 ]; then
    if [ -z "${SLURM_JOB_ID:-}" ]; then
        echo "[REFUSED] --run only inside a SLURM job (RI-HFXk memory); see header."; exit 1
    fi
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
    mpiexec -n "$SLURM_NTASKS" -bind-to "core:$OMP_NUM_THREADS" cp2k.psmp -i rihfx_kp.inp -o rihfx_kp.out
    if grep -q "PROGRAM ENDED AT" rihfx_kp.out; then
        echo "  [run] RI-HFX + k-points: RAN ($(grep -m1 'ENERGY| Total FORCE_EVAL' rihfx_kp.out | awk '{print $NF}') Ha)"
    else
        echo "  [run] RI-HFX + k-points: FAILED -- $(grep -m1 -A2 ABORT rihfx_kp.out | tr -s ' ' | tail -1)"
    fi
fi
echo "== outputs kept in $WORK"
