#!/usr/bin/env python3
"""
parse_surface_energy.py -- steps 06/07: surface energies and Wulff shape.

  gamma = (E_slab - n_fu * E_bulk_per_fu) / (2 A)        [J/m^2]

E_bulk_per_fu comes from 03_cellopt/bulk.json (same CUTOFF/REL_CUTOFF/basis
and k-density as the slabs -- checked). n_fu from 06's slab_meta.json, A from
07's run_meta.json (the area after rescaling to the CP2K lattice). Only
slabs whose GEO_OPT printed OPTIMIZATION COMPLETED are used. For each plane
the lowest-gamma termination enters the Wulff construction (pymatgen
WulffShape on the CP2K bulk lattice), whose facet area fractions are the
morphology.

  python ~/work_cp2k/parsers/parse_surface_energy.py --slabs-dir slabs_studies \\
      --bulk-json 03_cellopt/bulk.json
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import cp2k_blocks as cb  # noqa: E402
from core import cp2k_output as co  # noqa: E402

EV_A2_TO_J_M2 = 16.02176634


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slabs-dir", required=True)
    p.add_argument("--bulk-json", required=True)
    p.add_argument("--no-wulff", action="store_true")
    a = p.parse_args()

    bulk = json.loads(Path(a.bulk_json).read_text())
    rows = []
    for d in sorted(Path(a.slabs_dir).glob("slab_*")):
        meta_f, run_f = d / "slab_meta.json", d / "run_meta.json"
        if not (meta_f.exists() and run_f.exists()):
            continue
        meta, run = json.loads(meta_f.read_text()), json.loads(run_f.read_text())
        out = co.find_output(d)
        t = co.read(out) if out else ""
        status = "converged" if co.opt_converged(t) else ("running/failed" if t else "not run")
        for key in ("cutoff", "rel_cutoff", "basis"):
            if run.get(key) != bulk.get(key):
                print(f"[WARNING] {d.name}: {key} {run.get(key)} != bulk {bulk.get(key)}")
        row = dict(slab=d.name, hkl=meta["hkl"], termination=meta["termination"], status=status,
                   natoms=meta["natoms"], n_fu=meta["n_formula_units"], area=run["area_A2"])
        if status == "converged":
            e = co.final_energy_ev(t)
            row["E_eV"] = e
            row["gamma_J_m2"] = (e - meta["n_formula_units"] * bulk["E_per_fu_eV"]) / (2 * run["area_A2"]) * EV_A2_TO_J_M2
        rows.append(row)
    print(f"{'slab':<16} {'status':<15} {'gamma [J/m2]':>12}")
    for r in rows:
        g = f"{r['gamma_J_m2']:.4f}" if "gamma_J_m2" in r else "-"
        print(f"{r['slab']:<16} {r['status']:<15} {g:>12}")
    best = {}
    for r in rows:
        if "gamma_J_m2" in r:
            k = tuple(r["hkl"])
            if k not in best or r["gamma_J_m2"] < best[k]["gamma_J_m2"]:
                best[k] = r
    result = dict(bulk=a.bulk_json, slabs=rows,
                  best_per_plane={str(k): v["gamma_J_m2"] for k, v in best.items()})
    if best and not a.no_wulff and len(best) >= 2:
        from pymatgen.analysis.wulff import WulffShape
        lat = cb.read_structure(bulk["relaxed_cif"]).lattice
        ws = WulffShape(lat, list(best), [v["gamma_J_m2"] for v in best.values()])
        frac = {str(k): v for k, v in ws.area_fraction_dict.items()}
        result.update(wulff_area_fraction=frac, wulff_weighted_gamma=ws.weighted_surface_energy,
                      wulff_shape_factor=ws.shape_factor)
        print("\nWulff area fractions:")
        for k, v in sorted(frac.items(), key=lambda kv: -kv[1]):
            print(f"  {k:<14} {v * 100:6.1f} %")
        try:
            import matplotlib
            matplotlib.use("Agg")
            ws.get_plot().figure.savefig(Path(a.slabs_dir) / "wulff.png", dpi=150)
        except Exception as exc:
            print(f"[INFO] Wulff plot skipped: {exc}")
    (Path(a.slabs_dir) / "surface_energies.json").write_text(json.dumps(result, indent=1, default=str))


if __name__ == "__main__":
    main()
