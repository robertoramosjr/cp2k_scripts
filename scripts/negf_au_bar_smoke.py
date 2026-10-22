#!/usr/bin/env python3
"""Write a minimal CP2K NEGF input for a perfect Au [001] bar (smoke test, no molecule).

Device = 7 identical cells along z (L2 L1 | L0 | S | R0 | R1 R2), each cell = --repeats x a
with a 2x2-atom cross-section (4 Au per a, a = 4.078 A). L1/L2 (R1/R2) are the two bulk unit
cells of a contact, L0 (R0) the screening cell, S the scattering region. For a perfect bar
T(E) must be a flat integer-like staircase (number of open channels), which validates the
contact/screening/scattering bookkeeping before any molecule is added.
Template: tests/QS/regtest-negf-1/H28.inp of the CP2K tree; manual &NEGF.

Usage: negf_au_bar_smoke.py --output-dir DIR [--repeats 2] [--vacuum 14] [--cutoff 300]
"""
import argparse
from pathlib import Path

A_AU = 4.078
LABELS = ["L2", "L1", "L0", "S", "R0", "R1", "R2"]


def coords(repeats):
    rows = []
    for ic, lab in enumerate(LABELS):
        for r in range(repeats):
            z0 = (ic * repeats + r) * A_AU
            for x, y, dz in ((0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)):
                rows.append((x * A_AU / 2, y * A_AU / 2, z0 + dz * A_AU / 2, lab))
    return rows


def idx(label, repeats):
    n = 4 * repeats
    s = LABELS.index(label) * n + 1
    return f"{s}..{s + n - 1}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", default=".")
    ap.add_argument("--project", default="AuBar")
    ap.add_argument("--repeats", type=int, default=2, help="a-repeats per cell (cell length = repeats*a)")
    ap.add_argument("--vacuum", type=float, default=14.0, help="lateral box (A)")
    ap.add_argument("--cutoff", type=float, default=300.0)
    ap.add_argument("--rel-cutoff", type=float, default=40.0)
    ap.add_argument("--basis", default="DZVP-MOLOPT-SR-GTH")
    ap.add_argument("--temperature", type=float, default=300.0)
    a = ap.parse_args()

    rows = coords(a.repeats)
    lz = len(LABELS) * a.repeats * A_AU
    cd = "\n".join(f"      Au {x:12.6f} {y:12.6f} {z:12.6f} {lab}" for x, y, z, lab in rows)
    inp = f"""&GLOBAL
  PROJECT {a.project}
  RUN_TYPE NEGF
  PRINT_LEVEL LOW
  EXTENDED_FFT_LENGTHS .TRUE.
  PREFERRED_DIAG_LIBRARY SL
&END GLOBAL

&NEGF
  ENERGY_LBOUND -1.0
  EPS_DENSITY 1e-5
  INTEGRATION_METHOD SIMPSON
  MAX_SCF 100
  NPROC_POINT 1
  V_SHIFT_MAX_ITERS 100
  &CONTACT
    ELECTRIC_POTENTIAL [eV] 0.0
    REFINE_FERMI_LEVEL
    TEMPERATURE [K] {a.temperature}
    &BULK_REGION
      LIST {idx('L2', a.repeats)} {idx('L1', a.repeats)}
      &CELL
        LIST {idx('L2', a.repeats)}
      &END CELL
      &CELL
        LIST {idx('L1', a.repeats)}
      &END CELL
    &END BULK_REGION
    &SCREENING_REGION
      LIST {idx('L0', a.repeats)}
    &END SCREENING_REGION
  &END CONTACT
  &CONTACT
    ELECTRIC_POTENTIAL [eV] 0.0
    REFINE_FERMI_LEVEL
    TEMPERATURE [K] {a.temperature}
    &BULK_REGION
      LIST {idx('R1', a.repeats)} {idx('R2', a.repeats)}
      &CELL
        LIST {idx('R1', a.repeats)}
      &END CELL
      &CELL
        LIST {idx('R2', a.repeats)}
      &END CELL
    &END BULK_REGION
    &SCREENING_REGION
      LIST {idx('R0', a.repeats)}
    &END SCREENING_REGION
  &END CONTACT
  &MIXING
    ALPHA 0.2
    METHOD BROYDEN_MIXING
    NBUFFER 8
  &END MIXING
  &PRINT
    &DOS
      FILENAME device
      FROM_ENERGY [eV] -3
      TILL_ENERGY [eV] 3
      N_GRIDPOINTS 121
    &END DOS
    &PROGRAM_RUN_INFO LOW
    &END PROGRAM_RUN_INFO
    &TRANSMISSION
      FILENAME transm
      FROM_ENERGY [eV] -3
      TILL_ENERGY [eV] 3
      N_GRIDPOINTS 121
    &END TRANSMISSION
  &END PRINT
  &SCATTERING_REGION
    LIST {idx('S', a.repeats)}
  &END SCATTERING_REGION
&END NEGF

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF {a.cutoff}
      REL_CUTOFF {a.rel_cutoff}
      NGRIDS 5
    &END MGRID
    &QS
      EPS_DEFAULT 1E-12
    &END QS
    &POISSON
      PERIODIC XYZ
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1E-6
      MAX_SCF 100
      ADDED_MOS 10
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.3
        NBUFFER 8
      &END MIXING
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] {a.temperature}
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
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {a.vacuum} {a.vacuum} {lz:.6f}
      PERIODIC XYZ
    &END CELL
    &COORD
{cd}
    &END COORD
    &KIND Au
      BASIS_SET {a.basis}
      POTENTIAL GTH-PBE-q11
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{a.project}.inp").write_text(inp)
    print(f"wrote {out / (a.project + '.inp')}: {len(rows)} Au, cell {a.vacuum} x {a.vacuum} x {lz:.2f} A")


if __name__ == "__main__":
    main()
