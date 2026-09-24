#!/usr/bin/env python3
"""
03_cellopt.py -- production CELL_OPT, single step (no separate GEO_OPT).

CELL_OPT relaxes cell and ions together, so its final geometry and energy
are the production bulk (and the bulk reference for surface energies).
&CELL_REF (default 1.15x the lattice vectors) keeps the grids built on a
reference cell larger than any cell the optimizer visits; STRESS_TENSOR
ANALYTICAL; BFGS; MAX_FORCE in eV/angstrom; &SCF%PRINT%RESTART ON.

  python ~/work_cp2k/03_cellopt.py --structure bulk.vasp --basis DZVP-MOLOPT-SR-GTH \\
      --project ZrO2_m_cellopt --cutoff 700 --rel-cutoff 60 --kpoints 6 6 6
  sbatch ~/work_cp2k/jobs/job_single.sh $PWD
  python ~/work_cp2k/parsers/parse_cellopt.py --run-dir .   # -> relaxed CIF + bulk.json
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import cp2k_blocks as cb  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cb.add_common_args(p)
    p.add_argument("--cell-ref-factor", type=float, default=cb.DEFAULTS["cell_ref_factor"])
    p.add_argument("--max-iter", type=int, default=cb.DEFAULTS["max_iter"])
    p.add_argument("--max-force", type=float, default=cb.DEFAULTS["max_force"])
    p.add_argument("--pressure-tolerance", type=float, default=cb.DEFAULTS["pressure_tolerance"])
    p.add_argument("--keep-space-group", action="store_true",
                   help="MOTION/CELL_OPT KEEP_SPACE_GROUP (preserve the symmetry during the run).")
    a = p.parse_args()

    st = cb.read_structure(a.structure)
    text = cb.build_input(
        st, a, run_type="CELL_OPT", kpoints=a.kpoints, cell_ref_factor=a.cell_ref_factor,
        stress=True, scf_kwargs=dict(restart_print=True),
        motion_kwargs=dict(max_iter=a.max_iter, max_force=a.max_force,
                           pressure_tolerance=a.pressure_tolerance,
                           keep_space_group=a.keep_space_group))
    path = cb.write_input(text, a.output_dir, a.project)
    (Path(a.output_dir) / "run_meta.json").write_text(json.dumps(dict(
        step="03_cellopt", kpoints=a.kpoints, cutoff=a.cutoff, rel_cutoff=a.rel_cutoff,
        basis=a.basis, natoms=len(st), formula=st.composition.formula,
        k_density=cb.kmesh_density(st.lattice, a.kpoints), source=str(Path(a.structure).resolve())),
        indent=1))
    print(f"{st.composition.formula} ({len(st)} atoms) CELL_OPT -> {path}")


if __name__ == "__main__":
    main()
