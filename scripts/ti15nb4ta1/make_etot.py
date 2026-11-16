#!/usr/bin/env python3
"""Write the production total-energy (Etot) single point for the Ti15Nb4Ta1 bulk.

RUN_TYPE ENERGY_FORCE (ENERGY alone does not compute the stress) on the final geometry of a converged CELL_OPT (read from its -1.restart), with the chosen
k-mesh, full SCF (MAX_SCF 1500, no IGNORE_CONVERGENCE_FAILURE), STRESS_TENSOR ANALYTICAL printed (residual
pressure check), and NO &CELL_REF: the PW grid is built on the relaxed cell itself (CELL_REF lesson).
All other DFT settings come from bulk_pipeline.make_input (single source of blocks).

Usage (from Ti15Nb4Ta1/):
  python3 make_etot.py --restart bulk/03_kmesh/k6-5-5/TiNbTa_k6-5-5-1.restart --kmesh 6 5 5 \
      --outdir bulk/05_etot_k6-5-5
"""
import argparse
import json
import re
from pathlib import Path

import bulk_pipeline as bp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--restart", required=True, help="CELL_OPT -1.restart with the relaxed cell/coords")
    ap.add_argument("--kmesh", type=int, nargs=3, required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--decisions", default=str(bp.HERE / "bulk" / "decisions.json"), help="CUTOFF/REL_CUTOFF source")
    ap.add_argument("--project", default="TiNbTa_etot")
    a = ap.parse_args()

    dec = json.loads(Path(a.decisions).read_text())
    lat, sp, frac = bp.parse_restart(a.restart)
    # CELL_OPT template (full SCF + STRESS_TENSOR ANALYTICAL), then turned into a single point
    txt = bp.make_input(a.project, "CELL_OPT", lat, sp, frac, dec["cutoff"], dec["rel_cutoff"], kmesh=a.kmesh)
    txt = txt.replace("RUN_TYPE CELL_OPT", "RUN_TYPE ENERGY_FORCE")  # ENERGY alone skips the stress tensor
    txt = txt[:txt.index("&MOTION")]
    txt = txt.replace("&END SUBSYS\n&END FORCE_EVAL",
                      "&END SUBSYS\n  &PRINT\n    &STRESS_TENSOR\n    &END STRESS_TENSOR\n  &END PRINT\n&END FORCE_EVAL")
    assert "CELL_REF" not in txt and "MAX_SCF 1500" in txt and "IGNORE_CONVERGENCE" not in txt
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{a.project}.inp").write_text(txt)
    meta = dict(source_restart=str(Path(a.restart).resolve()), kmesh=a.kmesh, cutoff=dec["cutoff"],
                rel_cutoff=dec["rel_cutoff"], natoms=len(sp), volume_A3=abs(bp.det3(lat)),
                abc=[bp.norm(v) for v in lat])
    (out / "etot_meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
