#!/usr/bin/env python3
"""Write a CP2K PBE-D3(BJ) input for an Au tip / azobenzene / Au tip junction.

Reads an xyz written by build_azo_au_bridge.py (cell "a b c" in the comment line, Au tips plus
the molecule) and writes <outdir>/<project>.inp for GEO_OPT or ENERGY. Gamma point, OT-free
diagonalisation with Fermi-Dirac smearing (Au is metallic), PBE + D3(BJ), DZVP-MOLOPT-SR-GTH.
Only the --free-au Au atoms closest to the molecule on each side relax; the rest of the
electrodes stay fixed at the nanowire geometry. CUTOFF/REL_CUTOFF are required on purpose:
they must come from the grid convergence of this system (cp2k-pipeline, step 01), not a default.

Usage:
  make_junction_input.py --xyz structures/trans_phys.xyz --project trans_phys \
      --outdir trans_phys --cutoff 500 --rel-cutoff 50 [--run-type GEO_OPT] [--free-au 4]
  grid scan point (ENERGY, one SCF step, convergence failure ignored): add --scf-probe
"""
import argparse
from pathlib import Path

import numpy as np

BASIS = {"Au": "DZVP-MOLOPT-SR-GTH", "C": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH",
         "N": "DZVP-MOLOPT-SR-GTH", "S": "DZVP-MOLOPT-SR-GTH"}
POTENTIAL = {"Au": "GTH-PBE-q11", "C": "GTH-PBE-q4", "H": "GTH-PBE-q1",
             "N": "GTH-PBE-q5", "S": "GTH-PBE-q6"}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    cell = [float(x) for x in lines[1].split()[1:4]]
    atoms = [(l.split()[0], np.array([float(x) for x in l.split()[1:4]])) for l in lines[2:2 + n]]
    return cell, atoms


def free_au_indices(atoms, n_free):
    """1-based indices of the n_free Au atoms nearest to the molecule on each side (by z)."""
    mol = np.array([p for e, p in atoms if e != "Au"])
    zc = mol[:, 2].mean()
    free = []
    for side in (-1, 1):
        cand = [(np.min(np.linalg.norm(mol - p, axis=1)), i + 1)
                for i, (e, p) in enumerate(atoms) if e == "Au" and side * (p[2] - zc) > 0]
        free += [i for _, i in sorted(cand)[:n_free]]
    return sorted(free)


def fmt_list(idx):
    return " ".join(str(i) for i in idx)


def build(args, cell, atoms):
    elements = sorted({e for e, _ in atoms})
    unknown = [e for e in elements if e not in BASIS]
    if unknown:
        raise SystemExit(f"no basis/potential defined for {unknown}")
    au = [i + 1 for i, (e, _) in enumerate(atoms) if e == "Au"]
    free = set(free_au_indices(atoms, args.free_au))
    fixed = [i for i in au if i not in free]
    coord = "\n".join(f"      {e:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for e, p in atoms)
    kinds = "\n".join(
        f"    &KIND {e}\n      BASIS_SET {BASIS[e]}\n      POTENTIAL {POTENTIAL[e]}\n    &END KIND"
        for e in elements)
    motion = ""
    if args.run_type == "GEO_OPT":
        motion = f"""
&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER {args.max_iter}
    MAX_FORCE 4.5E-4
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      LIST {fmt_list(fixed)}
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
"""
    ignore = "\n      IGNORE_CONVERGENCE_FAILURE TRUE" if args.scf_probe else ""
    return f"""&GLOBAL
  PROJECT {args.project}
  RUN_TYPE {args.run_type}
  PRINT_LEVEL MEDIUM
  EXTENDED_FFT_LENGTHS .TRUE.
  PREFERRED_DIAG_LIBRARY SL
&END GLOBAL
{motion}
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF {args.cutoff}
      REL_CUTOFF {args.rel_cutoff}
      NGRIDS 5
    &END MGRID
    &QS
      EPS_DEFAULT 1E-12
      EXTRAPOLATION USE_GUESS
    &END QS
    &POISSON
      PERIODIC XYZ
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1E-6{ignore}
      MAX_SCF {1 if args.scf_probe else 100}
      ADDED_MOS 10
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.2
        NBUFFER 8
      &END MIXING
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300.0
      &END SMEAR
      &PRINT
        &RESTART
          FILENAME RESTART
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {cell[0]:.6f} {cell[1]:.6f} {cell[2]:.6f}
      PERIODIC XYZ
    &END CELL
    &COORD
{coord}
    &END COORD
{kinds}
  &END SUBSYS
&END FORCE_EVAL
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--project", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cutoff", type=float, required=True)
    ap.add_argument("--rel-cutoff", type=float, required=True)
    ap.add_argument("--run-type", choices=["GEO_OPT", "ENERGY"], default="GEO_OPT")
    ap.add_argument("--free-au", type=int, default=4, help="relaxing Au atoms per tip")
    ap.add_argument("--scf-probe", action="store_true",
                    help="MAX_SCF 1 + IGNORE_CONVERGENCE_FAILURE (grid convergence scan only)")
    ap.add_argument("--max-iter", type=int, default=300)
    args = ap.parse_args()
    cell, atoms = read_xyz(args.xyz)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{args.project}.inp"
    path.write_text(build(args, cell, atoms))
    n_au = sum(e == "Au" for e, _ in atoms)
    print(f"wrote {path}: {len(atoms)} atoms ({n_au} Au), cell {cell}")


if __name__ == "__main__":
    main()
