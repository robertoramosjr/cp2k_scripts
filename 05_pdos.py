#!/usr/bin/env python3
"""
05_pdos.py -- projected DOS on a supercell, Gamma-only + OT.

A Gamma-only supercell samples the Brillouin zone of the primitive cell by
folding, which lets the cheap OT minimizer do the SCF. &PRINT%PDOS with
COMPONENTS (per-l/m channels) and NLUMO -1 (all virtual states). Broadening
is done afterwards in Python (parsers/parse_pdos.py --npoints --sigma), not
with an input keyword.

--dos-interface:
  legacy   &DFT%PRINT%PDOS  (present in the CP2K 2026.1 schema on this cluster)
  unified  &DFT%PRINT%DOS with &PDOS nested (the 2026.2 layout seen on coaraci);
           run scripts/probe_build.sh first -- 2026.1 rejects it.

  python ~/work_cp2k/05_pdos.py --structure ZrO2_m_relaxed.cif --basis DZVP-MOLOPT-SR-GTH \\
      --project ZrO2_m_pdos --cutoff 700 --rel-cutoff 60 --supercell 3 3 3
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import cp2k_blocks as cb  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cb.add_common_args(p, need_kpoints=False)
    p.add_argument("--supercell", type=int, nargs=3, required=True, metavar=("NX", "NY", "NZ"))
    p.add_argument("--nlumo", type=int, default=-1)
    p.add_argument("--dos-interface", choices=["legacy", "unified"], default="legacy")
    a = p.parse_args()

    st = cb.read_structure(a.structure)
    sc = st.copy()
    sc.make_supercell(a.supercell)
    def pdos(ind):
        return (f"{ind}&PDOS\n{ind}  COMPONENTS\n{ind}  NLUMO {a.nlumo}\n{ind}&END PDOS\n")

    if a.dos_interface == "legacy":
        dft_print = "    &PRINT\n" + pdos(" " * 6) + "    &END PRINT\n"
    else:
        dft_print = "    &PRINT\n      &DOS\n" + pdos(" " * 8) + "      &END DOS\n    &END PRINT\n"
    text = cb.build_input(sc, a, run_type="ENERGY", kpoints=[1, 1, 1], dft_print=dft_print,
                          scf_kwargs=dict(restart_print=False))
    path = cb.write_input(text, a.output_dir, a.project)
    (Path(a.output_dir) / "run_meta.json").write_text(json.dumps(dict(
        step="05_pdos", supercell=a.supercell, natoms=len(sc), nlumo=a.nlumo,
        dos_interface=a.dos_interface), indent=1))
    print(f"supercell {a.supercell} -> {len(sc)} atoms, Gamma-only OT -> {path}")


if __name__ == "__main__":
    main()
