#!/usr/bin/env python3
"""
parse_pdos.py -- step 05: Gaussian-broadened PDOS from CP2K *.pdos files.

Each <project>-k<N>-1.pdos holds, per atomic kind, the eigenvalues [a.u.],
occupations and projections on each (l, m) component (COMPONENTS). The
header carries E(Fermi) [a.u.]. Broadening is done here (like VASP's
NEDOS=6000 grid + smearing), not with an input keyword.

Writes pdos.csv (energy relative to E_F, per kind and per l channel, and the
total) and pdos.png.

  python ~/work_cp2k/parsers/parse_pdos.py --run-dir 05_pdos --npoints 6000 --sigma 0.1
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.cp2k_blocks import HARTREE_EV  # noqa: E402

RE_HEAD = re.compile(r"atomic kind\s+(\S+).*E\(Fermi\)\s*=\s*([-\d.Ee+]+)\s*a\.u\.")


def read_pdos(path: Path):
    lines = path.read_text().splitlines()
    m = RE_HEAD.search(lines[0])
    kind, ef = m.group(1), float(m.group(2)) * HARTREE_EV
    cols = lines[1].lstrip("#").replace("Eigenvalue [a.u.]", "Eigenvalue").split()
    cols = cols[cols.index("Occupation") + 1:]
    data = np.loadtxt(path, comments="#")
    return kind, ef, data[:, 1] * HARTREE_EV, cols, data[:, 3:]


def channel(c: str) -> str:
    return c[0] if c[0] in "spdf" else c


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", default=".")
    p.add_argument("--npoints", type=int, default=6000)
    p.add_argument("--sigma", type=float, default=0.1, help="Gaussian width [eV].")
    p.add_argument("--emin", type=float, default=-10.0, help="Plot window relative to E_F [eV].")
    p.add_argument("--emax", type=float, default=10.0)
    a = p.parse_args()

    d = Path(a.run_dir)
    files = sorted(d.glob("*.pdos"))
    if not files:
        sys.exit("ERROR: no *.pdos in --run-dir.")
    parsed = [read_pdos(f) for f in files]
    ef = parsed[0][1]
    allE = np.concatenate([x[2] for x in parsed]) - ef
    grid = np.linspace(allE.min() - 5 * a.sigma, allE.max() + 5 * a.sigma, a.npoints)
    norm = 1.0 / (a.sigma * np.sqrt(2 * np.pi))

    def broaden(eps, w):
        return (w[None, :] * np.exp(-0.5 * ((grid[:, None] - eps[None, :]) / a.sigma) ** 2)).sum(1) * norm

    curves = {}
    for kind, _, eps, cols, proj in parsed:
        for l in sorted({channel(c) for c in cols}):
            idx = [i for i, c in enumerate(cols) if channel(c) == l]
            curves[f"{kind}_{l}"] = broaden(eps - ef, proj[:, idx].sum(1))
    curves["total"] = sum(v for k, v in curves.items())
    with open(d / "pdos.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E-EF_eV"] + list(curves))
        for i, e in enumerate(grid):
            w.writerow([f"{e:.4f}"] + [f"{curves[k][i]:.6f}" for k in curves])
    print(f"E_F = {ef:.4f} eV; {len(curves) - 1} channels -> {d / 'pdos.csv'}")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        for k, v in curves.items():
            ax.plot(grid, v, lw=1.6 if k == "total" else 0.9, label=k, c="k" if k == "total" else None)
        ax.axvline(0, c="gray", ls="--", lw=0.6)
        ax.set_xlim(a.emin, a.emax)
        ax.set_xlabel("E - E_F [eV]")
        ax.set_ylabel("PDOS [states/eV]")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(d / "pdos.png", dpi=150)
    except Exception as exc:
        print(f"[INFO] plot skipped: {exc}")


if __name__ == "__main__":
    main()
