#!/usr/bin/env python3
"""
04_bands.py -- PBE band structure: RUN_TYPE ENERGY + &DFT%PRINT%BAND_STRUCTURE.

Path from pymatgen HighSymmKpath (Setyawan-Curtarolo). The special points are
fractional coordinates of the reciprocal lattice of kpath.prim, so the input
cell IS kpath.prim (CP2K reads SPECIAL_POINT in UNITS B_VECTOR of the input
cell). The SCF mesh is step 3's mesh; if kpath.prim differs from the step-3
cell, the mesh is converted at constant density n_i*|a_i| (printed).

  python ~/work_cp2k/04_bands.py --structure ZrO2_m_relaxed.cif --basis DZVP-MOLOPT-SR-GTH \\
      --project ZrO2_m_bands --cutoff 700 --rel-cutoff 60 --kpoints 6 6 6
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import cp2k_blocks as cb  # noqa: E402

from pymatgen.symmetry.bandstructure import HighSymmKpath  # noqa: E402


def label(s: str) -> str:
    s = s.replace("\\Gamma", "GAMMA").replace("\\Sigma", "SIGMA").replace("\\", "")
    return re.sub(r"[^A-Za-z0-9_]", "", s) or "K"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cb.add_common_args(p)
    p.add_argument("--npoints", type=int, default=20, help="k-points between consecutive special points.")
    p.add_argument("--symprec", type=float, default=0.01)
    a = p.parse_args()

    st = cb.read_structure(a.structure)
    kp = HighSymmKpath(st, symprec=a.symprec)
    prim = kp.prim
    same = np.allclose(sorted(prim.lattice.abc), sorted(st.lattice.abc), atol=1e-3)
    mesh = list(a.kpoints) if same else cb.kmesh_from_density(
        prim.lattice, cb.kmesh_density(st.lattice, a.kpoints))
    if not same:
        print(f"[INFO] kpath.prim ({len(prim)} atoms) differs from the input cell ({len(st)} atoms); "
              f"SCF mesh {a.kpoints} -> {mesh} at constant density.")
    pts = kp.kpath["kpoints"]
    sets = []
    for branch in kp.kpath["path"]:
        lines = ["        &KPOINT_SET", "          UNITS B_VECTOR"]
        for name in branch:
            x, y, z = pts[name]
            lines.append(f"          SPECIAL_POINT {label(name)} {x:.8f} {y:.8f} {z:.8f}")
        lines += [f"          NPOINTS {a.npoints}", "        &END KPOINT_SET"]
        sets.append("\n".join(lines))
    dft_print = ("    &PRINT\n      &BAND_STRUCTURE\n"
                 f"        FILE_NAME {a.project}.bs\n        ADDED_MOS {a.added_mos}\n"
                 + "\n".join(sets) + "\n      &END BAND_STRUCTURE\n    &END PRINT\n")
    text = cb.build_input(prim, a, run_type="ENERGY", kpoints=mesh, dft_print=dft_print,
                          scf_kwargs=dict(restart_print=False))
    path = cb.write_input(text, a.output_dir, a.project)
    (Path(a.output_dir) / "run_meta.json").write_text(json.dumps(dict(
        step="04_bands", path=kp.kpath["path"], scf_kpoints=mesh, npoints=a.npoints,
        natoms=len(prim)), indent=1))
    print(f"path: {' | '.join('-'.join(b) for b in kp.kpath['path'])}\n-> {path}")


if __name__ == "__main__":
    main()
