#!/usr/bin/env python3
"""
07_slab_opt.py -- GEO_OPT input for a slab cut by 06, on the CP2K bulk lattice.

The slab from 06 was cut from the VASP bulk. Before relaxing it in CP2K, its
in-plane vectors are rebuilt from the SAME bulk translations [uvw] (stored in
slab_meta.json) on the CP2K-relaxed lattice (--bulk-cp2k, from
parsers/parse_cellopt.py), and heights along the normal are scaled by
d_hkl(CP2K)/d_hkl(ref). Without this the slab would carry the VASP->CP2K
lattice mismatch as artificial strain, biasing gamma.

Cell fixed (GEO_OPT), c along z, SURFACE_DIPOLE_CORRECTION along Z (zero for
the symmetric slabs 06 writes, kept as a safeguard), k-mesh
ceil(k_density/|a|) x ceil(k_density/|b|) x 1 with k_density = n_i*|a_i| of
the converged bulk mesh (step 03 run_meta.json: "k_density").

  python ~/work_cp2k/07_slab_opt.py --slab-dir slabs_studies/slab_-111_1 \\
      --bulk-cp2k 03_cellopt/ZrO2_m_relaxed.cif --basis DZVP-MOLOPT-SR-GTH \\
      --cutoff 700 --rel-cutoff 60 --k-density 30
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import cp2k_blocks as cb  # noqa: E402

from pymatgen.core import Lattice, Structure  # noqa: E402


def rescale(slab: Structure, meta: dict, bulk_ref: Structure, bulk_new: Structure):
    u1, u2 = (np.array(u) for u in meta["uvw_in_plane"])
    n1, n2 = u1 @ bulk_new.lattice.matrix, u2 @ bulk_new.lattice.matrix
    hkl = meta["hkl"]
    zscale = bulk_new.lattice.d_hkl(hkl) / bulk_ref.lattice.d_hkl(hkl)
    l1, l2 = np.linalg.norm(n1), np.linalg.norm(n2)
    gam = math.degrees(math.acos(np.dot(n1, n2) / (l1 * l2)))
    z = slab.cart_coords[:, 2]
    zc = 0.5 * (z.max() + z.min())
    t_new = (z.max() - z.min()) * zscale
    c = t_new + meta["vacuum_A"]
    lat = Lattice.from_parameters(l1, l2, c, 90.0, 90.0, gam)
    frac = slab.frac_coords.copy()
    frac[:, 2] = ((z - zc) * zscale + c / 2.0) / c
    strain = (l1 / slab.lattice.a - 1.0, l2 / slab.lattice.b - 1.0)
    return Structure(lat, slab.species, frac), strain, zscale, t_new


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slab-dir", required=True, help="Folder from 06 (POSCAR + slab_meta.json).")
    p.add_argument("--bulk-cp2k", default=None,
                   help="CP2K-relaxed bulk (same setting as the cut bulk). Omit ONLY for previews.")
    p.add_argument("--k-density", type=float, required=True)
    p.add_argument("--no-dipole", action="store_true")
    p.add_argument("--max-iter", type=int, default=cb.DEFAULTS["max_iter"])
    p.add_argument("--max-force", type=float, default=cb.DEFAULTS["max_force"])
    cb.add_common_args(p, need_kpoints=False)
    # structure/project/output-dir default from the slab folder
    for act in p._actions:
        if act.dest in ("structure", "project"):
            act.required = False
    p.set_defaults(mixing_alpha=0.2)  # vacuum -> charge sloshing; gentler than bulk 0.4
    a = p.parse_args()

    d = Path(a.slab_dir)
    meta = json.loads((d / "slab_meta.json").read_text())
    slab = cb.read_structure(a.structure or d / "POSCAR")
    a.project = a.project or f"slab_{''.join(str(h) for h in meta['hkl'])}_{meta['termination']}"
    a.output_dir = a.output_dir if a.output_dir != "." else str(d)
    if a.bulk_cp2k:
        if None in meta["uvw_in_plane"]:
            sys.exit("ERROR: slab_meta.json has no rational [uvw]; cannot rescale.")
        slab, strain, zscale, t_new = rescale(slab, meta, cb.read_structure(meta["bulk"]),
                                              cb.read_structure(a.bulk_cp2k))
    else:
        strain, zscale, t_new = (0.0, 0.0), 1.0, meta["thickness_A"]
        print("[WARNING] no --bulk-cp2k: slab left on the VASP lattice (preview only).")
    kx, ky, _ = cb.kmesh_from_density(slab.lattice, a.k_density, periodic=(True, True, False))
    text = cb.build_input(
        slab, a, run_type="GEO_OPT", kpoints=[kx, ky, 1], surface_dipole=not a.no_dipole,
        scf_kwargs=dict(restart_print=True),
        motion_kwargs=dict(max_iter=a.max_iter, max_force=a.max_force))
    path = cb.write_input(text, a.output_dir, a.project)
    slab.to(filename=str(Path(a.output_dir) / "slab_cp2k_start.cif"))
    run_meta = dict(step="07_slab_opt", kpoints=[kx, ky, 1], k_density=a.k_density,
                    cutoff=a.cutoff, rel_cutoff=a.rel_cutoff, basis=a.basis,
                    bulk_cp2k=str(Path(a.bulk_cp2k).resolve()) if a.bulk_cp2k else None,
                    inplane_strain_vs_vasp=list(strain), zscale=zscale,
                    area_A2=float(np.linalg.norm(np.cross(*slab.lattice.matrix[:2]))),
                    thickness_A=t_new, natoms=len(slab))
    (Path(a.output_dir) / "run_meta.json").write_text(json.dumps(run_meta, indent=1))
    print(f"{a.project}: {slab.composition.formula} ({len(slab)} at) k={kx}x{ky}x1 "
          f"strain={strain[0]*100:+.2f}%/{strain[1]*100:+.2f}% -> {path}")


if __name__ == "__main__":
    main()
