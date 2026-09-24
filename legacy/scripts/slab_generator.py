#!/usr/bin/env python3
"""
slab_generator.py

Step 6 of the pipeline (after Etot): turns an already-cut slab (e.g. a legacy
VASP POSCAR from Materials Studio/pymatgen) into a CP2K GEO_OPT input that is
consistent with the CP2K bulk reference, so surface energies

    gamma = (E_slab - n_fu * E_bulk_per_fu) / (2 * A)

compare slab and bulk on the SAME potential-energy surface.

What it does, in order:
  1. Finds the vacuum axis (largest empty gap) and refuses 2D-confined
     objects (wires/rods: vacuum along two axes) unless --allow-wire.
  2. Identifies the Miller plane by matching the two in-plane slab vectors to
     integer combinations [uvw] of --bulk-ref (the bulk the slab was cut from);
     hkl = [uvw]_1 x [uvw]_2. Rotation-invariant (only lengths + angle).
  3. If --bulk-cp2k is given (CIF from extract_relaxed_structure.py on the
     Etot restart), rebuilds the in-plane vectors from the SAME [uvw] on the
     CP2K lattice and scales the out-of-plane spacing by d_hkl(cp2k)/d_hkl(ref)
     -- removes the VASP->CP2K lattice mismatch as artificial strain.
  4. Reorients the cell: a in x, b in the xy plane, c exactly along z with
     length thickness + --vacuum. SURFACE_DIPOLE_CORRECTION (SURF_DIP_DIR Z)
     is only meaningful with the surface normal along z; legacy cells here
     had tilted normals (in-plane vectors with z components).
  5. k-mesh: kx, ky = ceil(k_density / |a|), ceil(k_density / |b|), kz = 1,
     where k_density = n_i * |a_i| from the converged BULK mesh.
  6. Writes <project>.inp (via input_generator.render_input, GEO_OPT, fixed
     cell), slab_input.cif and slab_meta.json (plane, n formula units, area,
     stoichiometry...) for the surface-energy parser.

Usage (one slab):
    python ~/work_cp2k/scripts/slab_generator.py \\
        --slab legacy/.../slab_1/POSCAR --bulk-ref legacy/.../CONTCAR \\
        --bulk-cp2k ZrO2_monoclinic_relaxed.cif \\
        --cutoff 500 --rel-cutoff 60 --k-density 30 \\
        --project-name ZrO2_m_001_1 --output-dir slabs_studies/slab_001_1
"""

import argparse
import itertools
import json
import math
import sys
from functools import reduce
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import input_generator as ig  # noqa: E402

from pymatgen.core import Lattice, Structure  # noqa: E402

UVW_MAX = 4
LEN_TOL = 0.03    # angstrom, in-plane vector length match against bulk-ref
ANGLE_TOL = 0.5   # degrees
WIRE_GAP = 5.0    # angstrom; a second empty gap this large means not a slab


def parse_args():
    p = argparse.ArgumentParser(description="Prepare a CP2K GEO_OPT slab input from a cut slab.")
    p.add_argument("--slab", required=True, help="Slab structure file (POSCAR, CIF...).")
    p.add_argument("--bulk-ref", required=True,
                   help="Bulk the slab was cut from (used to identify [uvw]/hkl).")
    p.add_argument("--bulk-cp2k", default=None,
                   help="CP2K-relaxed bulk (same phase/setting as --bulk-ref). If given, the "
                        "slab is rescaled onto this lattice. Omit only for previews/cost estimates.")
    p.add_argument("--vacuum", type=float, default=15.0, help="Vacuum thickness along z (angstrom).")
    p.add_argument("--k-density", type=float, required=True,
                   help="n_i*|a_i| (angstrom) of the converged bulk mesh; kx,ky = ceil(k_density/|a|).")
    p.add_argument("--allow-wire", action="store_true",
                   help="Do not abort when vacuum is found along two axes.")
    p.add_argument("--no-surface-dipole", action="store_true",
                   help="Skip SURFACE_DIPOLE_CORRECTION (only sensible for symmetric slabs).")
    p.add_argument("--project-name", required=True)
    p.add_argument("--output-dir", required=True)
    # Passed through to input_generator.render_input -- same defaults as there.
    p.add_argument("--functional", default="PBE", choices=["PBE", "PBE0", "HSE06"])
    p.add_argument("--basis-family", default=ig.BASIS_FAMILY_DEFAULT)
    p.add_argument("--cutoff", type=int, required=True)
    p.add_argument("--rel-cutoff", type=int, default=60)
    p.add_argument("--basis-set-file", default="BASIS_MOLOPT")
    p.add_argument("--potential-file", default="GTH_POTENTIALS")
    p.add_argument("--kind-overrides", default=None)
    p.add_argument("--eps-scf", type=float, default=1.0E-6)
    p.add_argument("--max-scf", type=int, default=300)
    p.add_argument("--outer-max-scf", type=int, default=50)
    p.add_argument("--mixing-alpha", type=float, default=0.2,
                   help="Lower than the bulk default (0.4): slabs with vacuum are prone to "
                        "charge sloshing.")
    p.add_argument("--max-iter", type=int, default=500)
    p.add_argument("--max-force", type=float, default=0.010)
    return p.parse_args()


def vacuum_gaps(s):
    """Largest empty gap (angstrom) along each lattice direction, and where it ends."""
    rec = s.lattice.reciprocal_lattice_crystallographic.matrix
    out = []
    for i in range(3):
        d = 1.0 / np.linalg.norm(rec[i])
        f = np.sort(s.frac_coords[:, i] % 1.0)
        gaps = np.diff(np.concatenate([f, [f[0] + 1.0]]))
        k = int(np.argmax(gaps))
        out.append((gaps[k] * d, (f[k] + gaps[k]) % 1.0))
    return out


def gcd_reduce(v):
    g = reduce(math.gcd, [abs(int(x)) for x in v if x != 0], 0)
    return tuple(int(x // g) for x in v) if g else tuple(int(x) for x in v)


def match_uvw(bulk, t1, t2):
    """All ([uvw]1, [uvw]2) bulk lattice-vector pairs reproducing |t1|, |t2| and their angle."""
    A = bulk.lattice.matrix
    l1, l2 = np.linalg.norm(t1), np.linalg.norm(t2)
    ang = math.degrees(math.acos(np.dot(t1, t2) / (l1 * l2)))
    cands = [(np.array(u), np.dot(u, A)) for u in itertools.product(range(-UVW_MAX, UVW_MAX + 1), repeat=3)
             if u != (0, 0, 0)]
    c1 = [(u, v) for u, v in cands if abs(np.linalg.norm(v) - l1) < LEN_TOL]
    c2 = [(u, v) for u, v in cands if abs(np.linalg.norm(v) - l2) < LEN_TOL]
    pairs = []
    for u1, v1 in c1:
        for u2, v2 in c2:
            a = math.degrees(math.acos(np.clip(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1)))
            if abs(a - ang) < ANGLE_TOL and np.linalg.norm(np.cross(u1, u2)) > 0:
                pairs.append((u1, u2))
    return pairs


def d_spacing(lattice, hkl):
    return lattice.d_hkl(hkl)


def main():
    args = parse_args()
    slab = Structure.from_file(args.slab)
    bulk_ref = Structure.from_file(args.bulk_ref)
    bulk_cp2k = Structure.from_file(args.bulk_cp2k) if args.bulk_cp2k else None

    gaps = vacuum_gaps(slab)
    iv = int(np.argmax([g[0] for g in gaps]))
    others = [i for i in range(3) if i != iv]
    extra = [gaps[i][0] for i in others if gaps[i][0] > WIRE_GAP]
    if extra and not args.allow_wire:
        sys.exit(f"ERROR: vacuum along two axes ({gaps[iv][0]:.1f} and {extra[0]:.1f} A) -- "
                 f"this is a wire/rod, not a slab. Use --allow-wire to force.")

    L = slab.lattice.matrix
    t1, t2 = L[others[0]], L[others[1]]
    if np.dot(np.cross(t1, t2), L[iv]) < 0:  # keep a right-handed (t1, t2, normal) frame
        t1, t2 = t2, t1
        others = others[::-1]
    pairs = match_uvw(bulk_ref, t1, t2)
    if not pairs:
        sys.exit("ERROR: in-plane vectors do not match any [uvw] pair of --bulk-ref "
                 "(wrong bulk, or slab already strained?).")
    u1, u2 = pairs[0]
    hkl = gcd_reduce(np.cross(u1, u2))

    # Unwrap across the vacuum so the slab is contiguous, then get Cartesian coords.
    fc = slab.frac_coords.copy()
    fc[:, iv] = (fc[:, iv] - gaps[iv][1] + 1e-4) % 1.0
    cart = slab.lattice.get_cartesian_coords(fc)

    # Orthonormal frame: x along t1, z along the surface normal.
    ex = t1 / np.linalg.norm(t1)
    ez = np.cross(t1, t2); ez /= np.linalg.norm(ez)
    ey = np.cross(ez, ex)
    R = np.vstack([ex, ey, ez])
    xyz = cart @ R.T
    # In-plane fractional coords w.r.t. (t1, t2), height along the normal.
    inplane = np.linalg.solve(np.array([[t1 @ ex, t2 @ ex], [t1 @ ey, t2 @ ey]]), xyz[:, :2].T).T
    z = xyz[:, 2] - xyz[:, 2].min()

    if bulk_cp2k is not None:
        n1, n2 = np.dot(u1, bulk_cp2k.lattice.matrix), np.dot(u2, bulk_cp2k.lattice.matrix)
        zscale = d_spacing(bulk_cp2k.lattice, hkl) / d_spacing(bulk_ref.lattice, hkl)
        strain = (np.linalg.norm(n1) / np.linalg.norm(t1) - 1, np.linalg.norm(n2) / np.linalg.norm(t2) - 1)
    else:
        n1, n2, zscale, strain = t1, t2, 1.0, (0.0, 0.0)
    z = z * zscale

    l1, l2 = np.linalg.norm(n1), np.linalg.norm(n2)
    gam = math.acos(np.dot(n1, n2) / (l1 * l2))
    thickness = float(z.max())
    cz = thickness + args.vacuum
    a_new = [l1, 0.0, 0.0]
    b_new = [l2 * math.cos(gam), l2 * math.sin(gam), 0.0]
    c_new = [0.0, 0.0, cz]
    lat = Lattice([a_new, b_new, c_new])
    frac = np.column_stack([inplane[:, 0], inplane[:, 1], (z + args.vacuum / 2) / cz])
    new = Structure(lat, slab.species, frac)

    kx = max(1, math.ceil(args.k_density / l1))
    ky = max(1, math.ceil(args.k_density / l2))

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    run_args = argparse.Namespace(
        project_name=args.project_name, run_type="GEO_OPT", functional=args.functional,
        basis_family=args.basis_family, cutoff=args.cutoff, rel_cutoff=args.rel_cutoff,
        stress_rel_cutoff_factor=1.0, eps_scf=args.eps_scf, max_scf=args.max_scf,
        outer_max_scf=args.outer_max_scf, mixing_alpha=args.mixing_alpha,
        max_iter=args.max_iter, max_force=args.max_force, kpoints=[kx, ky, 1],
        basis_set_file=args.basis_set_file, potential_file=args.potential_file,
        kind_overrides=args.kind_overrides, surface_dipole=not args.no_surface_dipole,
    )
    (out / f"{args.project_name}.inp").write_text(ig.render_input(new, run_args), encoding="utf-8")
    new.to(filename=str(out / "slab_input.cif"))

    comp = new.composition
    bulk_comp = bulk_ref.composition.reduced_composition
    n_fu = comp.get_reduced_composition_and_factor()[1]
    stoich = comp.reduced_composition == bulk_comp
    meta = dict(
        source=str(Path(args.slab).resolve()), bulk_ref=str(Path(args.bulk_ref).resolve()),
        bulk_cp2k=str(Path(args.bulk_cp2k).resolve()) if args.bulk_cp2k else None,
        hkl=list(hkl), uvw_in_plane=[list(map(int, u1)), list(map(int, u2))],
        formula=comp.formula, natoms=len(new), stoichiometric=bool(stoich),
        n_formula_units=float(n_fu) if stoich else None,
        area_A2=float(np.linalg.norm(np.cross(a_new, b_new))),
        thickness_A=thickness, vacuum_A=args.vacuum, cell_c_A=cz,
        inplane_strain_vs_ref=[float(x) for x in strain], kpoints=[kx, ky, 1],
        cutoff=args.cutoff, rel_cutoff=args.rel_cutoff,
        surface_dipole_correction=not args.no_surface_dipole,
    )
    (out / "slab_meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")

    print(f"{args.project_name}: hkl={hkl} {comp.formula} ({len(new)} atoms), "
          f"A={meta['area_A2']:.2f} A^2, thickness={thickness:.2f} A, c={cz:.2f} A, "
          f"k={kx}x{ky}x1, strain={strain[0]*100:+.2f}%/{strain[1]*100:+.2f}%"
          + ("" if stoich else "  [NON-STOICHIOMETRIC: gamma needs mu_O]"))


if __name__ == "__main__":
    main()
