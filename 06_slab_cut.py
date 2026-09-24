#!/usr/bin/env python3
"""
06_slab_cut.py -- cut symmetric, stoichiometric, non-polar slabs from a bulk.

For surface energies gamma = (E_slab - n*E_bulk_fu) / (2A) to mean the energy
of ONE surface, both faces must be the same termination (symmetric slab),
the slab must have the bulk stoichiometry (else gamma depends on mu_O), and
it must carry no net dipole along the normal (Tasker type I/II). This script
enumerates every termination pymatgen's SlabGenerator finds for each plane
(plain, symmetrized, and Tasker III -> II reconstructions of polar ones),
keeps only candidates that pass all three tests, removes duplicates, and
ranks them by broken cation-anion bonds per surface area (lowest first, a
standard proxy for the most stable termination). Candidate n is written to
<output-dir>/slab_<hkl>_<n>/ with c orthogonal to the surface and along z
(what SURFACE_DIPOLE_CORRECTION / SURF_DIP_DIR Z in 07 needs).

Miller indices refer to the lattice of --bulk as given (e.g. the standard
P2_1/c setting, beta > 90, for m-ZrO2: (-111) and (111) are different planes).

  python ~/work_cp2k/06_slab_cut.py --bulk CONTCAR --planes="-1,1,1;1,1,1" \\
      --min-slab 12 --vacuum 15 --output-dir slabs_studies
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import cp2k_blocks as cb  # noqa: E402

from pymatgen.analysis.structure_matcher import StructureMatcher  # noqa: E402
from pymatgen.core import Lattice, Structure  # noqa: E402
from pymatgen.core.surface import SlabGenerator  # noqa: E402


def hkl_tag(hkl) -> str:
    return "".join(str(int(h)) for h in hkl)


def parse_planes(s: str):
    return [tuple(int(x) for x in t.split(",")) for t in s.split(";") if t.strip()]


def count_bonds(st: Structure, cation: str, anion: str, rcut: float) -> int:
    n = 0
    for i, site in enumerate(st):
        if site.specie.element.symbol != cation:
            continue
        n += sum(1 for nb in st.get_neighbors(site, rcut) if nb.specie.element.symbol == anion)
    return n


def thickness(slab) -> float:
    n = np.cross(slab.lattice.matrix[0], slab.lattice.matrix[1])
    n /= np.linalg.norm(n)
    z = slab.cart_coords @ n
    return float(z.max() - z.min())


def to_z_oriented(slab, vacuum: float) -> Structure:
    """c orthogonal to the surface, a along x, c along z; slab centred, vacuum reset."""
    s = slab.get_orthogonal_c_slab()
    L = s.lattice
    a, b = L.matrix[0], L.matrix[1]
    n = np.cross(a, b)
    n /= np.linalg.norm(n)
    z = s.cart_coords @ n
    t = z.max() - z.min()
    c = t + vacuum
    lat = Lattice.from_parameters(L.a, L.b, c, 90.0, 90.0, L.gamma)
    frac = s.frac_coords.copy()
    frac[:, 2] = (z - z.min() + vacuum / 2.0) / c
    return Structure(lat, s.species, frac)


def tasker3_to_2(slab, layer_tol: float = 0.3, max_combos: int = 5000):
    """Tasker III -> II by hand: move half of one outermost layer to the bulk
    position one stacking period beyond the opposite face.

    A stoichiometric polar slab made of n stacked oriented unit cells (OUC)
    has atoms at r + j*t (j = 0..n-1, t = OUC stacking vector). Moving an atom
    of the top layer by -n*t puts it exactly where the next bulk layer below
    the bottom face would be (and vice versa), so the result keeps the bulk
    stoichiometry and bulk-like sites. Every half-subset of the outermost
    layer is tried on the 1x1 cell and on 2x1/1x2/2x2 supercells (a layer of
    one atom cannot be halved on 1x1); returns the non-polar symmetric ones.
    """
    from itertools import combinations

    from pymatgen.core.surface import Slab

    ouc = slab.oriented_unit_cell
    n_rep = len(slab) / len(ouc)
    if abs(n_rep - round(n_rep)) > 1e-6:
        return []
    out = []
    for sc in ([1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1]):
        s = slab.copy()
        if sc != [1, 1, 1]:
            s.make_supercell(sc)
        t = ouc.lattice.matrix[2]
        nrm = np.cross(s.lattice.matrix[0], s.lattice.matrix[1])
        nrm /= np.linalg.norm(nrm)
        z = s.cart_coords @ nrm
        shift = round(n_rep) * t * (1 if np.dot(t, nrm) > 0 else -1)
        for side, sign in (("top", -1), ("bottom", +1)):
            edge = z.max() if side == "top" else z.min()
            layer = [i for i in range(len(s)) if abs(z[i] - edge) < layer_tol]
            species = {s[i].specie.element.symbol for i in layer}
            if len(species) != 1 or len(layer) % 2:
                continue
            combos = combinations(layer, len(layer) // 2)
            for k, move in enumerate(combos):
                if k >= max_combos:
                    break
                cart = s.cart_coords.copy()
                for i in move:
                    cart[i] = cart[i] + sign * shift
                new = Slab(s.lattice, s.species, cart, s.miller_index, s.oriented_unit_cell,
                           s.shift, s.scale_factor, coords_are_cartesian=True,
                           reorient_lattice=False)
                if not new.is_polar() and new.is_symmetric():
                    out.append(new)
        if out:
            break
    return out


def uvw_of(bulk: Structure, vec, max_den: int = 4) -> list[float] | None:
    """[uvw] with vec = u*a + v*b + w*c of the bulk (same Cartesian frame; the
    generator runs with reorient_lattice=False for this). Rational, not only
    integer: a primitive fcc translation is (1/2, 1/2, 0) of the conventional
    cell. 07 rebuilds the in-plane vectors as uvw @ lattice(CP2K bulk)."""
    x = np.linalg.solve(bulk.lattice.matrix.T, vec)
    for den in range(1, max_den + 1):
        r = np.round(x * den) / den
        if np.allclose(x, r, atol=1e-4):
            return [float(v) for v in r]
    return None


def hkl_from_normal(bulk: Structure, normal) -> tuple[int, ...]:
    """Miller indices of the plane with Cartesian normal `normal`: h_i = a_i . n,
    reduced to the smallest integers (self-check of the cut)."""
    from fractions import Fraction
    from functools import reduce
    from math import gcd
    h = bulk.lattice.matrix @ normal
    h = h / np.abs(h[np.abs(h) > 1e-6]).min()
    fr = [Fraction(x).limit_denominator(12) for x in h]
    den = reduce(lambda x, y: x * y // gcd(x, y), [f.denominator for f in fr])
    ints = [int(f * den) for f in fr]
    g = reduce(gcd, [abs(i) for i in ints if i])
    return tuple(i // g for i in ints)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bulk", required=True, help="Converged bulk (e.g. VASP CONTCAR), setting defines hkl.")
    p.add_argument("--planes", required=True,
                   help='Miller indices; write --planes="-1,1,1;1,1,1" (the "=" is needed '
                        'when the first index is negative).')
    p.add_argument("--min-slab", type=float, default=12.0, help="Minimum slab thickness (angstrom).")
    p.add_argument("--vacuum", type=float, default=15.0, help="Vacuum along z (angstrom).")
    p.add_argument("--max-terminations", type=int, default=3,
                   help="Write at most this many valid terminations per plane (ranked).")
    p.add_argument("--cation", default="Zr")
    p.add_argument("--anion", default="O")
    p.add_argument("--bond-cutoff", type=float, default=2.7, help="Cation-anion bond cutoff (angstrom).")
    p.add_argument("--tol", type=float, default=0.1, help="SlabGenerator layer tolerance (angstrom).")
    p.add_argument("--oxidation", default="Zr:4,O:-2",
                   help="Formal charges for the polarity test. Without them pymatgen's "
                        "Slab.is_polar() is always False.")
    p.add_argument("--output-dir", default=".")
    a = p.parse_args()

    bulk = cb.read_structure(a.bulk)
    bulk.add_oxidation_state_by_element(
        {k: float(v) for k, v in (x.split(":") for x in a.oxidation.split(","))})
    bulk_red = bulk.composition.reduced_composition
    nfu_bulk = cb.n_formula_units(bulk)
    bonds_per_fu = count_bonds(bulk, a.cation, a.anion, a.bond_cutoff) / nfu_bulk
    print(f"bulk {bulk.composition.formula}: {bonds_per_fu:.2f} {a.cation}-{a.anion} bonds per f.u.")
    matcher = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=2, primitive_cell=False, scale=False)
    summary = []

    for hkl in parse_planes(a.planes):
        gen = SlabGenerator(bulk, hkl, a.min_slab, a.vacuum, lll_reduce=True, center_slab=True,
                            primitive=True, max_normal_search=max(abs(h) for h in hkl) + 1,
                            reorient_lattice=False)  # keep the bulk Cartesian frame (uvw_of)
        raw = gen.get_slabs(tol=a.tol) + gen.get_slabs(tol=a.tol, symmetrize=True)
        cands = []
        for s in raw:
            cands.append(("plain", s))
            if s.is_polar():
                try:
                    cands += [("tasker2", t) for t in s.get_tasker2_slabs()]
                except Exception:
                    pass
                if s.composition.reduced_composition == bulk_red:
                    cands += [("tasker3to2", t) for t in tasker3_to_2(s)]
        valid = []
        for how, s in cands:
            stoich = s.composition.reduced_composition == bulk_red
            if not (stoich and not s.is_polar() and s.is_symmetric()):
                continue
            if any(matcher.fit(s, v["slab"]) for v in valid):
                continue
            nfu = cb.n_formula_units(s)
            broken = nfu * bonds_per_fu - count_bonds(s, a.cation, a.anion, a.bond_cutoff)
            area = s.surface_area
            valid.append(dict(slab=s, how=how, nfu=nfu, broken_per_A2=broken / (2 * area),
                              area=area, thick=thickness(s)))
        valid.sort(key=lambda v: (round(v["broken_per_A2"], 4), len(v["slab"])))
        tag = hkl_tag(hkl)
        if not valid:
            print(f"({tag}): NO symmetric stoichiometric non-polar termination found "
                  f"({len(cands)} candidates) -- try a larger --min-slab or --tol.")
            summary.append(dict(hkl=list(hkl), n_valid=0))
            continue
        for n, v in enumerate(valid[: a.max_terminations], start=1):
            s = v["slab"]
            out = Path(a.output_dir) / f"slab_{tag}_{n}"
            out.mkdir(parents=True, exist_ok=True)
            st = to_z_oriented(s, a.vacuum)
            st.remove_oxidation_states()
            st.to(filename=str(out / "POSCAR"), fmt="poscar")
            st.to(filename=str(out / "slab.cif"))
            ortho = s.get_orthogonal_c_slab()
            nrm = np.cross(ortho.lattice.matrix[0], ortho.lattice.matrix[1])
            hkl_chk = hkl_from_normal(bulk, nrm / np.linalg.norm(nrm))
            if hkl_chk not in (tuple(hkl), tuple(-x for x in hkl)):
                print(f"[WARNING] ({tag}) #{n}: normal gives {hkl_chk}, not {hkl}")
            meta = dict(
                hkl=list(hkl), termination=n, origin=v["how"], formula=st.composition.formula,
                natoms=len(st), n_formula_units=v["nfu"], stoichiometric=True, symmetric=True,
                polar=False, area_A2=float(v["area"]), thickness_A=v["thick"], vacuum_A=a.vacuum,
                broken_bonds_per_A2=v["broken_per_A2"], bond_cutoff_A=a.bond_cutoff,
                hkl_from_normal=list(hkl_chk),
                uvw_in_plane=[uvw_of(bulk, ortho.lattice.matrix[0]), uvw_of(bulk, ortho.lattice.matrix[1])],
                bulk=str(Path(a.bulk).resolve()), bulk_abc=list(bulk.lattice.abc),
                bulk_angles=list(bulk.lattice.angles), n_valid_terminations=len(valid))
            (out / "slab_meta.json").write_text(json.dumps(meta, indent=1))
            print(f"({tag}) #{n}: {st.composition.formula:<10} {len(st):>3} at, A={v['area']:6.2f} A2, "
                  f"t={v['thick']:5.2f} A, broken={v['broken_per_A2']:.4f}/A2 [{v['how']}] -> {out}")
        summary.append(dict(hkl=list(hkl), n_valid=len(valid),
                            written=min(len(valid), a.max_terminations)))
    (Path(a.output_dir) / "cut_summary.json").write_text(json.dumps(dict(
        bulk=str(Path(a.bulk).resolve()), min_slab=a.min_slab, vacuum=a.vacuum,
        planes=summary), indent=1))


if __name__ == "__main__":
    main()
