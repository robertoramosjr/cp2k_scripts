#!/usr/bin/env python3
"""
parse_growth_shape.py -- step 09: crystallite shape vs growth plane (Winterbottom).

Reads the gammas of parse_surface_energy.py (best termination per plane) and,
for each plane (hkl) taken as the growth / contact plane, builds the
Winterbottom shape: the Wulff polyhedron n_i . x <= gamma_i with the contact
facet moved to

  h_int = gamma_int - gamma_sub = gamma_hkl - W_adh = gamma_hkl (1 - beta)

beta = W_adh / gamma_hkl: 0 = no adhesion (full Wulff shape sitting on that
face), 1 = cut through the Wulff centre (half crystal), -> 2 = wetting film.
For each case it reports the height H along the growth normal, the
equivalent in-plane diameter D (projected area), H/D, the contact-area
fraction and the exposed facet fractions. H/D is scale free; H and D are
for V = 1000 A^3. Each shape is also written as OBJ/MTL in
<slabs-dir>/growth_shapes_obj/ (one polygon per facet, one material per
{hkl} family, contact face grey; Y-up, so Blender's default import stands the
crystal on its contact face) for rendering outside matplotlib.

  python ~/work_cp2k/cp2k_scripts/parsers/parse_growth_shape.py \\
      --slabs-dir slabs_studies --bulk-json 03_cellopt/bulk.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, HalfspaceIntersection

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import cp2k_blocks as cb  # noqa: E402

V_REF = 1000.0  # A^3


def shape(normals, offsets):
    """Convex polyhedron n_i . x <= d_i -> (vertices, hull)."""
    A = np.asarray(normals)
    d = np.asarray(offsets)
    # Chebyshev centre: max r s.t. n_i . x + r <= d_i
    res = linprog([0, 0, 0, -1], A_ub=np.hstack([A, np.ones((len(A), 1))]), b_ub=d,
                  bounds=[(None, None)] * 3 + [(0, None)])
    hs = HalfspaceIntersection(np.hstack([A, -d[:, None]]), res.x[:3])
    v = hs.intersections
    return v, ConvexHull(v)


def facet_areas(hull, normals, offsets, tol=1e-6):
    """Area per halfspace index (hull triangles matched by their plane)."""
    area = np.zeros(len(normals))
    for simp, eq in zip(hull.simplices, hull.equations):
        p = hull.points[simp]
        a = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        i = np.argmax(np.asarray(normals) @ eq[:3])
        if abs(np.dot(normals[i], p[0]) - offsets[i]) < 1e-4 * (1 + abs(offsets[i])):
            area[i] += a
    return area


def frame(n):
    """Orthonormal (e1, e2, n) with n the growth normal (pointing up)."""
    t = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(n, t)
    e1 /= np.linalg.norm(e1)
    return np.array([e1, np.cross(n, e1), n])


def write_obj(path, verts, normals, fam, contact, colors):
    """One Winterbottom shape as OBJ + MTL: one polygon per facet, one material
    per {hkl} family (+ "contact"). verts are in the growth frame (z up, contact
    on z = 0); written Y-up (OBJ convention), so Blender's default import stands
    the crystal on its contact face."""
    v = np.asarray(verts, float)
    hull = ConvexHull(v)
    faces = {}
    for simp, eq in zip(hull.simplices, hull.equations):
        faces.setdefault(int(np.argmax(normals @ eq[:3])), set()).update(int(j) for j in simp)
    used = sorted(set().union(*faces.values()))
    vid = {j: k + 1 for k, j in enumerate(used)}
    mtl = path.with_suffix(".mtl")
    names = {}
    lines = [f"mtllib {mtl.name}", f"o {path.stem}"]
    lines += [f"v {v[j, 0]:.6f} {v[j, 2]:.6f} {-v[j, 1]:.6f}" for j in used]
    for i in sorted(faces, key=lambda i: (fam[i], i)):
        name = "contact" if i == contact else "f" + "".join(str(x).replace("-", "m") for x in fam[i])
        names[name] = (0.6, 0.6, 0.6) if i == contact else colors[fam[i]][:3]
        ids = sorted(faces[i])
        c = v[ids].mean(axis=0)
        e1, e2, n = frame(normals[i])
        ang = [np.arctan2((v[j] - c) @ e2, (v[j] - c) @ e1) for j in ids]
        ids = [j for _, j in sorted(zip(ang, ids))]  # counter-clockwise seen from outside
        lines += [f"g {name}", f"usemtl {name}", "f " + " ".join(str(vid[j]) for j in ids)]
    path.write_text("\n".join(lines) + "\n")
    mtl.write_text("".join(f"newmtl {k}\nKd {r:.4f} {g:.4f} {b:.4f}\nKa 0 0 0\nKs 0.1 0.1 0.1\nd 1.0\n\n"
                           for k, (r, g, b) in names.items()))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slabs-dir", required=True)
    p.add_argument("--bulk-json", required=True)
    p.add_argument("--betas", type=float, nargs="+", default=[0.0, 0.5, 1.0, 1.5])
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--no-obj", action="store_true",
                   help="Skip the per-shape OBJ/MTL export (growth_shapes_obj/, for Blender etc.).")
    a = p.parse_args()

    from pymatgen.analysis.wulff import WulffShape
    sd = Path(a.slabs_dir)
    se = json.loads((sd / "surface_energies.json").read_text())
    bulk = json.loads(Path(a.bulk_json).read_text())
    lat = cb.read_structure(bulk["relaxed_cif"]).lattice
    best = {tuple(int(x) for x in k.strip("()").split(",")): g for k, g in se["best_per_plane"].items()}
    ws = WulffShape(lat, list(best), list(best.values()))
    fac = [f for f in ws.facets]
    N = np.array([f.normal / np.linalg.norm(f.normal) for f in fac])
    G = np.array([f.e_surf for f in fac])
    fam = [tuple(f.miller) for f in fac]
    families = list(best)

    results = []
    for hkl in families:
        # contact facet: the member of the family whose normal is "most -z"-like;
        # growth direction = -n_contact
        idx = [i for i, f in enumerate(fam) if f == hkl]
        ic = min(idx, key=lambda i: (N[i] @ np.array([0.1, 0.2, 1.0])))
        up = -N[ic]
        for beta in a.betas:
            d = G.copy()
            d[ic] = G[ic] * (1 - beta)
            v, hull = shape(N, d)
            s = (V_REF / hull.volume) ** (1 / 3)
            v = v * s
            hull = ConvexHull(v)
            area = facet_areas(hull, N, d * s)
            R = frame(up)
            loc = v @ R.T
            H = loc[:, 2].max() - loc[:, 2].min()
            proj = ConvexHull(loc[:, :2]).volume  # 2D hull "volume" = area
            D = 2 * np.sqrt(proj / np.pi)
            free = area.copy()
            free[ic] = 0
            byfam = {str(f): float(sum(free[i] for i in range(len(fam)) if fam[i] == f) / free.sum())
                     for f in families}
            top = max((i for i in range(len(N)) if i != ic), key=lambda i: N[i] @ up)
            results.append(dict(growth_plane=str(hkl), beta=beta, gamma_hkl=float(G[ic]),
                                H_A=float(H), D_eq_A=float(D), H_over_D=float(H / D),
                                contact_fraction=float(area[ic] / area.sum()),
                                top_facet=str(fam[top]) if N[top] @ up > 0.999 and area[top] > 1e-6 * area.sum() else "edge/vertex",
                                exposed_area_fraction=byfam,
                                _verts=loc.tolist(), _contact=int(ic), _hkl=hkl))

    print(f"Winterbottom shapes, V = {V_REF:.0f} A^3 (beta = W_adh / gamma_hkl)\n")
    print(f"{'plane':<11}{'gamma':>7}{'beta':>6}{'H [A]':>8}{'D [A]':>8}{'H/D':>7}{'contact':>9}  top      dominant exposed")
    for r in results:
        dom = sorted(r["exposed_area_fraction"].items(), key=lambda kv: -kv[1])[:2]
        dom = ", ".join(f"{k} {v * 100:.0f}%" for k, v in dom)
        print(f"{r['growth_plane']:<11}{r['gamma_hkl']:7.3f}{r['beta']:6.2f}{r['H_A']:8.2f}{r['D_eq_A']:8.2f}"
              f"{r['H_over_D']:7.2f}{r['contact_fraction'] * 100:8.1f}%  {r['top_facet']:<11}{dom}")
    out = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    (sd / "growth_shapes.json").write_text(json.dumps(dict(V_ref_A3=V_REF, gammas=se["best_per_plane"],
                                                           shapes=out), indent=1))
    if a.no_plot and a.no_obj:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    cols = dict(zip(families, plt.cm.tab10.colors))
    if not a.no_obj:
        od = sd / "growth_shapes_obj"
        od.mkdir(exist_ok=True)
        for r in results:
            v = np.array(r["_verts"])
            v[:, 2] -= v[:, 2].min()
            Nl = N @ frame(-N[r["_contact"]]).T
            tag = "".join(str(x).replace("-", "m") for x in r["_hkl"])
            write_obj(od / f"growth_{tag}_beta{r['beta']:.2f}.obj", v, Nl, fam, r["_contact"], cols)
        print(f"[INFO] {len(results)} OBJ/MTL in {od} (Y-up; facets grouped per family, contact = grey)")
    if a.no_plot:
        return
    nb = len(a.betas)
    elev, azim = np.radians(12), np.radians(35)
    view = np.array([np.cos(elev) * np.cos(azim), np.cos(elev) * np.sin(azim), np.sin(elev)])
    fig = plt.figure(figsize=(3.0 * nb, 3.0 * len(families)))
    for k, r in enumerate(results):
        ax = fig.add_subplot(len(families), nb, k + 1, projection="3d")
        v = np.array(r["_verts"])
        v[:, 2] -= v[:, 2].min()
        R = frame(-N[r["_contact"]])
        Nl = N @ R.T
        hull = ConvexHull(v)
        for simp, eq in zip(hull.simplices, hull.equations):
            if eq[:3] @ view <= 0:  # back face: mplot3d does not depth-sort
                continue
            i = np.argmax(Nl @ eq[:3])
            c = "0.6" if i == r["_contact"] else cols[fam[i]]
            ax.add_collection3d(Poly3DCollection([v[simp]], facecolor=c, edgecolor=c, lw=0.3, alpha=1.0))
        L = np.abs(v).max() * 1.05
        xx, yy = np.meshgrid([-L, L], [-L, L])
        ax.plot_surface(xx, yy, np.zeros_like(xx), color="0.85", alpha=0.4)
        ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(0, 2 * L)
        ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=np.degrees(elev), azim=np.degrees(azim)); ax.set_axis_off()
        ax.set_title(f"{r['growth_plane']}  β={r['beta']:.1f}\nH/D={r['H_over_D']:.2f}", fontsize=9)
    handles = [plt.Rectangle((0, 0), 1, 1, color=cols[f]) for f in families] + \
              [plt.Rectangle((0, 0), 1, 1, color="0.6")]
    fig.legend(handles, [f"{{{''.join(map(str, f))}}}  γ={best[f]:.3f}" for f in families] + ["contact"],
               loc="upper center", ncol=len(families) + 1, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(sd / "growth_shapes.png", dpi=130)
    print(f"\n[INFO] {sd / 'growth_shapes.png'}, {sd / 'growth_shapes.json'}")


if __name__ == "__main__":
    main()
