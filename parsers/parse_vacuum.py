#!/usr/bin/env python3
"""
parse_vacuum.py -- vacuum convergence of slab energies (single points from
07_slab_opt.py --run-type ENERGY --vacuum V, one folder per <slab>_vacVV).

For each slab, the surface-energy error of vacuum V is
    d_gamma(V) = [E(V) - E(V_max)] / (2A)          [J/m^2]
(bulk term cancels: same slab, same atoms). The recommended vacuum is the
smallest V such that |d_gamma| < --threshold for V and every larger V; the
global recommendation is the largest of the per-slab values.

  python ~/work_cp2k/parsers/parse_vacuum.py --scan-dir vacuum_convergence
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import cp2k_output as co  # noqa: E402

EV_A2_TO_J_M2 = 16.02176634
RE_DIR = re.compile(r"^(?P<slab>.+)_vac(?P<v>[\d.]+)$")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan-dir", required=True)
    p.add_argument("--threshold", type=float, default=0.005, help="J/m^2 (default 5 mJ/m^2).")
    a = p.parse_args()

    groups = defaultdict(list)
    for d in sorted(Path(a.scan_dir).iterdir()):
        m = RE_DIR.match(d.name)
        if not (d.is_dir() and m):
            continue
        out = co.find_output(d)
        t = co.read(out) if out else ""
        e = co.final_energy_ev(t) if co.ended(t) else None
        meta = json.loads((d / "run_meta.json").read_text())
        groups[m["slab"]].append(dict(v=float(m["v"]), E=e, area=meta["area_A2"], ranks=co.mpi_ranks(t)))
    summary, worst = {}, None
    for slab, rows in sorted(groups.items()):
        rows.sort(key=lambda r: r["v"])
        done = [r for r in rows if r["E"] is not None]
        if len(done) < 2 or done[-1]["v"] != rows[-1]["v"]:
            print(f"{slab}: incomplete ({len(done)}/{len(rows)} finished, need the largest vacuum)")
            continue
        ref = done[-1]
        for r in done:
            r["dgamma"] = (r["E"] - ref["E"]) / (2 * r["area"]) * EV_A2_TO_J_M2
        rec = next((r["v"] for i, r in enumerate(done)
                    if all(abs(x["dgamma"]) < a.threshold for x in done[i:])), None)
        summary[slab] = dict(recommended=rec, rows=done)
        print(f"{slab:<28} " + "  ".join(f"{r['v']:g}A:{r['dgamma'] * 1000:+.2f}" for r in done)
              + f"  (mJ/m2)  -> {rec} A")
        if rec is not None:
            worst = rec if worst is None else max(worst, rec)
    print(f"\nRecommended vacuum (all slabs, |d_gamma| < {a.threshold * 1000:g} mJ/m2): {worst} A")
    (Path(a.scan_dir) / "vacuum_summary.json").write_text(json.dumps(
        dict(threshold_J_m2=a.threshold, recommended=worst, slabs=summary), indent=1))


if __name__ == "__main__":
    main()
