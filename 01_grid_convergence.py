#!/usr/bin/env python3
"""
01_grid_convergence.py -- CUTOFF / REL_CUTOFF scan, official CP2K 2D method.

Follows the CP2K how-to "Converging the CUTOFF and REL_CUTOFF": each point is
RUN_TYPE ENERGY, Gamma-only, MAX_SCF 1 from the same ATOMIC guess, so the
energy differences between points isolate the grid error (the unconverged
SCF error is identical for every point and cancels). PRINT_LEVEL MEDIUM makes
CP2K print the "count for grid N" lines: how many Gaussian products were
mapped to each multigrid level. parsers/parse_grid.py reads both.

Two phases, run in order:
  --phase cutoff     vary CUTOFF at fixed REL_CUTOFF (--rel-cutoff, default 60)
  --phase relcutoff  vary REL_CUTOFF at the CUTOFF chosen in phase 1 (--cutoff)

Usage (from the phase folder, any cwd works):
  python ~/work_cp2k/01_grid_convergence.py --phase cutoff --structure bulk.vasp \\
      --basis DZVP-MOLOPT-SR-GTH --project ZrO2_c --values 300 400 500 600 700 800 900 1000
  python ~/work_cp2k/01_grid_convergence.py --phase relcutoff --cutoff 700 ... \\
      --values 30 40 50 60 70 80 100
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import cp2k_blocks as cb  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cb.add_common_args(p, need_kpoints=False)
    p.add_argument("--phase", required=True, choices=["cutoff", "relcutoff"])
    p.add_argument("--values", type=float, nargs="+", required=True,
                   help="CUTOFF values (phase cutoff) or REL_CUTOFF values (phase relcutoff), Ry.")
    a = p.parse_args()

    st = cb.read_structure(a.structure)
    root = Path(a.output_dir) / a.phase
    project = a.project
    for v in a.values:
        tag = f"{a.phase}{v:g}"
        if a.phase == "cutoff":
            a.cutoff = v
        else:
            a.rel_cutoff = v
        a.project = f"{project}_{tag}"
        text = cb.build_input(
            st, a, run_type="ENERGY", kpoints=[1, 1, 1], print_level="MEDIUM",
            scf_kwargs=dict(max_scf=1, outer_scf=False, restart_print=False,
                            ignore_convergence_failure=True))
        path = cb.write_input(text, root / tag, a.project)
        print(f"  {tag:<16} -> {path}")
    meta = dict(step="01_grid_convergence", phase=a.phase, values=a.values,
                natoms=len(st), formula=st.composition.formula, basis=a.basis,
                fixed_cutoff=a.cutoff if a.phase == "relcutoff" else None,
                fixed_rel_cutoff=a.rel_cutoff if a.phase == "cutoff" else None)
    (root / "scan_meta.json").write_text(json.dumps(meta, indent=1))
    print(f"\n{len(a.values)} inputs in {root}/  (atoms: {len(st)})")
    print(f"Next: sbatch --array=0-{len(a.values) - 1} ~/work_cp2k/jobs/job_scan_array.sh {root.resolve()}")
    print(f"Then: python ~/work_cp2k/parsers/parse_grid.py --scan-dir {root}")


if __name__ == "__main__":
    main()
