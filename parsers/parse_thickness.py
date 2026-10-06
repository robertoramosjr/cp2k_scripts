#!/usr/bin/env python3
"""
parse_thickness.py -- slab-thickness convergence of gamma (folders from
08_thickness_series.py: <slab>_thk<T>/, one plane/termination per <slab>).

For each thickness:
  gamma(N) = (E_slab(N) - N * E_bulk_per_fu) / (2A)        [J/m^2]
with E_bulk_per_fu from bulk.json (03). A bulk reference that is not exactly
consistent with the slabs makes gamma drift linearly with N, so the
Fiorentini-Methfessel fit E_slab(N) = 2A*gamma_FM + N*E_fit is also reported:
its slope E_fit vs bulk.json (meV/f.u.) is the consistency check, its intercept
gives gamma_FM (needs >= 3 thicknesses to mean anything).
The recommended thickness is the thinnest whose gamma, and that of every
thicker slab, is within --threshold of the thickest one.
Only GEO_OPT runs with OPTIMIZATION COMPLETED (or finished ENERGY runs) count.

  python ~/work_cp2k/parsers/parse_thickness.py --scan-dir thickness_convergence \\
      --bulk-json 03_cellopt/bulk.json
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import cp2k_output as co  # noqa: E402

EV_A2_TO_J_M2 = 16.02176634
RE_DIR = re.compile(r"^(?P<slab>.+)_thk(?P<t>[\d.]+)$")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan-dir", required=True)
    p.add_argument("--bulk-json", required=True)
    p.add_argument("--threshold", type=float, default=0.01, help="J/m^2 (default 10 mJ/m^2).")
    a = p.parse_args()

    bulk = json.loads(Path(a.bulk_json).read_text())
    groups = defaultdict(list)
    for d in sorted(Path(a.scan_dir).iterdir()):
        m = RE_DIR.match(d.name)
        if not (d.is_dir() and m and (d / "run_meta.json").exists()):
            continue
        run = json.loads((d / "run_meta.json").read_text())
        meta = json.loads((d / "slab_meta.json").read_text())
        for key in ("cutoff", "rel_cutoff", "basis"):
            if run.get(key) != bulk.get(key):
                print(f"[WARNING] {d.name}: {key} {run.get(key)} != bulk {bulk.get(key)}")
        out = co.find_output(d)
        t = co.read(out) if out else ""
        done = co.opt_converged(t) if run["run_type"] == "GEO_OPT" else (co.ended(t) and not co.scf_unconverged(t))
        row = dict(dir=d.name, thk_A=run["thickness_A"], natoms=run["natoms"], n_fu=meta["n_formula_units"],
                   area=run["area_A2"], status="ok" if done else ("running/failed" if t else "not run"))
        if done:
            row["E_eV"] = co.final_energy_ev(t)
            row["gamma_J_m2"] = (row["E_eV"] - row["n_fu"] * bulk["E_per_fu_eV"]) / (2 * row["area"]) * EV_A2_TO_J_M2
        groups[m["slab"]].append(row)

    result = {}
    for slab, rows in groups.items():
        rows.sort(key=lambda r: r["n_fu"])
        ok = [r for r in rows if "gamma_J_m2" in r]
        print(f"\n{slab}")
        print(f"  {'dir':<28} {'t [A]':>6} {'N_fu':>5} {'status':<15} {'gamma':>8} {'d_gamma':>8}")
        ref = ok[-1]["gamma_J_m2"] if ok else None
        for r in rows:
            g = r.get("gamma_J_m2")
            r["d_gamma_J_m2"] = None if g is None else g - ref
            gs = "-" if g is None else f"{g:8.4f}"
            ds = "-" if g is None else f"{g - ref:+8.4f}"
            print(f"  {r['dir']:<28} {r['thk_A']:6.2f} {r['n_fu']:5d} {r['status']:<15} {gs:>8} {ds:>8}")
        res = dict(points=rows)
        if len(ok) >= 2:
            rec = ok[-1]
            for r in reversed(ok):
                if all(abs(x["gamma_J_m2"] - ref) < a.threshold for x in ok if x["n_fu"] >= r["n_fu"]):
                    rec = r
            res["recommended"] = dict(dir=rec["dir"], thk_A=rec["thk_A"], n_fu=rec["n_fu"],
                                      note="relative to the thickest slab, which is only converged if the "
                                           "last points agree")
            print(f"  recommended (|d_gamma| < {a.threshold} J/m2): {rec['dir']} ({rec['thk_A']:.2f} A)")
        if len(ok) >= 3:
            n = np.array([r["n_fu"] for r in ok], float)
            e = np.array([r["E_eV"] for r in ok])
            slope, icpt = np.polyfit(n, e, 1)
            area = ok[-1]["area"]
            res["fiorentini_methfessel"] = dict(
                E_fit_per_fu_eV=slope, dE_vs_bulk_meV_per_fu=(slope - bulk["E_per_fu_eV"]) * 1e3,
                gamma_FM_J_m2=icpt / (2 * area) * EV_A2_TO_J_M2)
            print(f"  FM fit: E_fu {slope:.6f} eV ({(slope - bulk['E_per_fu_eV']) * 1e3:+.2f} meV/f.u. vs bulk.json), "
                  f"gamma_FM {icpt / (2 * area) * EV_A2_TO_J_M2:.4f} J/m2")
        result[slab] = res
    (Path(a.scan_dir) / "thickness_summary.json").write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
