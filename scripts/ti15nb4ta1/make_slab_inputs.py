#!/usr/bin/env python3
"""Write CP2K GEO_OPT inputs for the Ti15Nb4Ta1 slabs cut by make_slabs.py.

Consistent with the bulk reference (bulk/05_etot_k6-5-5): same DFT blocks (bulk_pipeline.make_input:
PBE, DZVP-MOLOPT-SR-GTH, CUTOFF/REL_CUTOFF, Gaussian smearing 0.01 eV, EPS_SCF, MAX_SCF 1500), plus
  * cell rotated so that a is along x, b in the xy plane and c (normal to the slab) along +z
    (SURFACE_DIPOLE_CORRECTION is implemented only for a normal parallel to a Cartesian axis);
  * SURFACE_DIPOLE_CORRECTION .TRUE., SURF_DIP_DIR Z (the disordered slab has inequivalent faces);
  * k-mesh ceil(dens/|a|) x ceil(dens/|b|) x 1, dens = max_i(n_i |a_i|) of the bulk 6-5-5 mesh;
  * Broyden alpha 0.2 (slab), ADDED_MOS 20 per 20-atom period;
  * GEO_OPT LBFGS, MAX_ITER 800, MAX_FORCE 0.01 eV/A (VASP EDIFFG), all atoms free;
  * two stages per slab: <proj>_s1.inp (ENERGY, no dipole, EPS_SCF 1E-5, writes the .kp wfn) then
    <proj>.inp (GEO_OPT, dipole at the fixed vacuum centre SURF_DIP_POS, SCF_GUESS RESTART,
    EXTRAPOLATION USE_PREV_P: USE_GUESS redoes the guess each geometry step and aborts the dipole).
    The default SURF_DIP_POS (min density, re-searched every step) jumps around in the flat vacuum and
    makes the SCF oscillate; a fixed position from the ATOMIC guess aborts at SCF step 2.
N per plane = smallest number of 20-atom periods with N*d >= --min-thk (default 9 A).

  python3 make_slab_inputs.py --slabs slabs_k655 --out slab_opt --cutoff 600 --rel 30
"""
import argparse, json, math, re
from pathlib import Path

import bulk_pipeline as bp

BULK_META = "bulk/05_etot_k6-5-5/etot_meta.json"


def read_poscar(path):
    L = Path(path).read_text().splitlines()
    sc = float(L[1].split()[0])
    lat = [[float(x) * sc for x in L[i].split()[:3]] for i in (2, 3, 4)]
    els, cnt = L[5].split(), [int(x) for x in L[6].split()]
    n = sum(cnt)
    frac = [[float(x) for x in L[8 + i].split()[:3]] for i in range(n)]
    sp = [e for e, c in zip(els, cnt) for _ in range(c)]
    return lat, sp, frac


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def unit(v):
    n = math.sqrt(dot(v, v))
    return [x / n for x in v]


def rotate(lat):
    """Rigid rotation: a -> x, b in xy, normal -> +z. Fractional coordinates are unchanged."""
    x = unit(lat[0])
    z = unit(cross(lat[0], lat[1]))
    y = cross(z, x)
    new = [[dot(v, x), dot(v, y), dot(v, z)] for v in lat]
    c = new[2]
    if abs(c[0]) > 1e-6 or abs(c[1]) > 1e-6 or c[2] <= 0:
        raise SystemExit(f"c is not along +z after rotation: {c}")
    return new


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slabs", default="slabs_k655")
    ap.add_argument("--out", default="slab_opt")
    ap.add_argument("--cutoff", type=float, default=600.0)
    ap.add_argument("--rel", type=float, default=30.0)
    ap.add_argument("--min-thk", type=float, default=9.0, help="minimum N*d [A]")
    ap.add_argument("--nplus", type=int, default=0, help="add this many periods (thickness series)")
    ap.add_argument("--term", default="t1")
    ap.add_argument("--big-natoms", type=int, default=60,
                    help="above this: 12 ranks x 4 threads + PARALLEL_GROUP_SIZE 0 (memory); else 24x2, -1")
    a = ap.parse_args()

    meta = json.loads(Path(BULK_META).read_text())
    dens = max(n * l for n, l in zip(meta["kmesh"], meta["abc"]))
    rows = []
    for hdir in sorted(Path(a.slabs).glob("hkl_*")):
        hkl = hdir.name[4:]
        cands = []
        for f in hdir.glob(f"slab_{hkl}_N*_{a.term}.vasp"):
            n = int(re.search(r"_N(\d+)_", f.name).group(1))
            cands.append((n, f))
        cands.sort()
        # period d = slope of the atomic span vs N (exact for slabs of the same termination)
        span = {}
        for n, f in cands:
            lat, sp, fr = read_poscar(f)
            cz = math.sqrt(dot(lat[2], lat[2]))
            z = sorted(r[2] * cz for r in fr)
            span[n] = z[-1] - z[0]
        ns = sorted(span)
        d = (span[ns[-1]] - span[ns[0]]) / (ns[-1] - ns[0])
        nsel = next((n for n in ns if n * d >= a.min_thk - 1e-6), None)
        if nsel is None:
            raise SystemExit(f"{hkl}: no N in {ns} reaches {a.min_thk} A (d={d:.2f}); cut more periods")
        nsel += a.nplus
        f = dict(cands)[nsel] if nsel in dict(cands) else None
        if f is None:
            raise SystemExit(f"{hkl}: N={nsel} not cut; rerun make_slabs.py --nperiods ... {nsel}")
        lat, sp, fr = read_poscar(f)
        if dot(cross(lat[0], lat[1]), lat[2]) < 0:   # left-handed cut: swap a<->b (same structure)
            lat = [lat[1], lat[0], lat[2]]
            fr = [[r[1], r[0], r[2]] for r in fr]
        lat = rotate(lat)
        la, lb = (math.sqrt(dot(v, v)) for v in lat[:2])
        kmesh = (math.ceil(dens / la), math.ceil(dens / lb), 1)
        proj = f"slab_{hkl}"
        nper = len(sp) // 20
        big = len(sp) > a.big_natoms
        # physical vacuum centre along z (SURF_DIP_POS is measured from z=0 of the cell, in A; verified by scan)
        cz = lat[2][2]
        zs = sorted((r[2] % 1.0) * cz for r in fr)
        zpos = ((zs[0] + zs[-1] + cz) / 2.0) % cz
        common = [
            ("      ADDED_MOS 20\n", f"      ADDED_MOS {20 * nper}\n"),
            ("        ALPHA 0.4\n", "        ALPHA 0.2\n"),
            ("      &PRINT\n        &RESTART OFF\n        &END RESTART\n      &END PRINT\n", ""),  # keep the .kp wfn
        ]
        if big:   # one k-point group over all ranks: matrices distributed (100-atom slab OOMs with -1)
            common.append(("      PARALLEL_GROUP_SIZE -1\n", "      PARALLEL_GROUP_SIZE 0\n"))

        def patch(txt, pairs):
            for old, new in pairs:
                assert txt.count(old) == 1, old
                txt = txt.replace(old, new)
            return txt
        # stage 1: no dipole, loose SCF from ATOMIC -> wavefunction (.kp). The first unmixed step from the
        # atomic guess spills density into the vacuum (>1E-5) and aborts the dipole correction otherwise.
        s1 = patch(bp.make_input(f"{proj}_s1", "GEO_OPT", lat, sp, fr, a.cutoff, a.rel, kmesh=kmesh),
                   common + [("  RUN_TYPE GEO_OPT\n", "  RUN_TYPE ENERGY\n"),
                             (f"      EPS_SCF {bp.EPS_SCF:.2E}\n", "      EPS_SCF 1.00E-05\n")])
        # stage 2: GEO_OPT with the dipole correction at the fixed vacuum centre, restarting the wfn
        txt = patch(bp.make_input(proj, "GEO_OPT", lat, sp, fr, a.cutoff, a.rel, kmesh=kmesh),
                    common + [("    &QS\n", "    SURFACE_DIPOLE_CORRECTION .TRUE.\n    SURF_DIP_DIR Z\n"
                               f"    SURF_DIP_POS {zpos:.3f}\n    &QS\n"),
                              ("    BASIS_SET_FILE_NAME BASIS_MOLOPT\n",
                               f"    BASIS_SET_FILE_NAME BASIS_MOLOPT\n    WFN_RESTART_FILE_NAME {proj}-RESTART.kp\n"),
                              ("      SCF_GUESS ATOMIC\n", "      SCF_GUESS RESTART\n"),
                              # USE_GUESS re-runs SCF_GUESS at every geometry step: the fresh guess spills
                              # density into the vacuum and aborts the dipole correction (job 854/855)
                              ("      EXTRAPOLATION USE_GUESS\n", "      EXTRAPOLATION USE_PREV_P\n")])
        txt += ("&MOTION\n  &GEO_OPT\n    OPTIMIZER LBFGS\n"
                f"    MAX_ITER {bp.MAX_ITER}\n    MAX_FORCE [eV*angstrom^-1] {bp.MAX_FORCE}\n  &END GEO_OPT\n"
                "  &PRINT\n    &RESTART\n      BACKUP_COPIES 3\n      &EACH\n        GEO_OPT 1\n      &END EACH\n"
                "    &END RESTART\n  &END PRINT\n&END MOTION\n")
        od = Path(a.out) / f"hkl_{hkl}_N{nsel}_{a.term}"
        od.mkdir(parents=True, exist_ok=True)
        (od / f"{proj}.inp").write_text(txt)
        (od / f"{proj}_s1.inp").write_text(s1)
        area = abs(cross(lat[0], lat[1])[2])
        info = dict(hkl=hkl, N=nsel, natoms=len(sp), d_A=d, thk_Nd_A=nsel * d, area_A2=area,
                    kmesh=list(kmesh), dens_A=dens, source=str(f), c_A=lat[2][2],
                    cutoff=a.cutoff, rel_cutoff=a.rel, bulk_ref=BULK_META,
                    surf_dip_pos_A=zpos, ranks_threads=[12, 4] if big else [24, 2],
                    parallel_group_size=0 if big else -1)
        (od / "slab_meta.json").write_text(json.dumps(info, indent=1))
        rows.append((od, proj, info))
        print(f"{hkl:>7} N={nsel} {len(sp):4d} at  N*d={nsel * d:5.2f} A  area={area:6.1f} A2  "
              f"k={kmesh[0]}x{kmesh[1]}x1  dip@{zpos:5.2f}/{cz:5.2f}  {'12x4 pg0' if big else '24x2'}  -> {od}")
    (Path(a.out) / "list_all.txt").write_text("".join(
        f"{od} {p} {i['ranks_threads'][0]} {i['ranks_threads'][1]}\n" for od, p, i in rows))


if __name__ == "__main__":
    main()
