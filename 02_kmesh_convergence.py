#!/usr/bin/env python3
"""
02_kmesh_convergence.py -- k-mesh convergence with a FULL CELL_OPT per mesh.

Each mesh relaxes cell + ions (VASP ISIF=3 analog) at the CUTOFF/REL_CUTOFF
from step 01, so the criterion is on the physical observables that matter:
relaxed energy per atom AND relaxed volume (parsers/parse_kmesh.py:
dE [meV/atom] and dV [%] against the densest mesh).

Meshes: --k-list "2,2,2;4,4,4" or --densities 20 30 40 (n_i*|a_i| in
angstrom, converted per lattice vector -- better for anisotropic cells).

  python ~/work_cp2k/02_kmesh_convergence.py --structure bulk.vasp \\
      --basis DZVP-MOLOPT-SR-GTH --project ZrO2_m --cutoff 700 --rel-cutoff 60 \\
      --densities 15 20 25 30 40
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
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--k-list", help='Explicit meshes, e.g. "2,2,2;4,4,4;6,6,6".')
    g.add_argument("--densities", type=float, nargs="+",
                   help="k-point densities n_i*|a_i| (angstrom); duplicates are dropped.")
    p.add_argument("--cell-ref-factor", type=float, default=cb.DEFAULTS["cell_ref_factor"])
    p.add_argument("--max-iter", type=int, default=cb.DEFAULTS["max_iter"])
    p.add_argument("--max-force", type=float, default=cb.DEFAULTS["max_force"])
    p.add_argument("--pressure-tolerance", type=float, default=cb.DEFAULTS["pressure_tolerance"])
    a = p.parse_args()

    st = cb.read_structure(a.structure)
    if a.k_list:
        meshes = [tuple(int(x) for x in t.split(",")) for t in a.k_list.split(";")]
    else:
        meshes = []
        for d in a.densities:
            m = tuple(cb.kmesh_from_density(st.lattice, d))
            if m not in meshes:
                meshes.append(m)
    root = Path(a.output_dir)
    project = a.project
    for m in meshes:
        tag = "k{}-{}-{}".format(*m)
        a.project = f"{project}_{tag}"
        text = cb.build_input(
            st, a, run_type="CELL_OPT", kpoints=list(m), cell_ref_factor=a.cell_ref_factor,
            stress=True, scf_kwargs=dict(restart_print=False),  # array: no concurrent .kp writes
            motion_kwargs=dict(max_iter=a.max_iter, max_force=a.max_force,
                               pressure_tolerance=a.pressure_tolerance))
        path = cb.write_input(text, root / tag, a.project)
        print(f"  {tag:<12} density={cb.kmesh_density(st.lattice, m):6.1f} A -> {path}")
    (root / "scan_meta.json").write_text(json.dumps(dict(
        step="02_kmesh_convergence", meshes=meshes, natoms=len(st), cutoff=a.cutoff,
        rel_cutoff=a.rel_cutoff, basis=a.basis, abc=st.lattice.abc), indent=1))
    print(f"\nNext: sbatch --array=0-{len(meshes) - 1} ~/work_cp2k/jobs/job_scan_array.sh {root.resolve()}")
    print(f"Then: python ~/work_cp2k/parsers/parse_kmesh.py --scan-dir {root}")


if __name__ == "__main__":
    main()
