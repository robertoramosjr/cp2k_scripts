"""
cp2k_output.py -- small, tolerant readers for CP2K text output shared by parsers/.

Convergence is never inferred from "PROGRAM ENDED AT" alone: CP2K does not
abort when the SCF or the optimizer fails, it prints the unconverged result.
"""

from __future__ import annotations

import re
from pathlib import Path

from .cp2k_blocks import HARTREE_EV

RE_ENERGY = re.compile(r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[(?:a\.u\.|hartree)\]:?\s+([-\d.Ee+]+)")
RE_GRID = re.compile(r"count for grid\s+(\d+):\s+(\d+)\s+cutoff \[a\.u\.\]\s+([-\d.Ee+]+)")
RE_NPROC = re.compile(r"Total number of message passing processes\s+(\d+)")
RE_VOLUME = re.compile(r"CELL\| Volume \[angstrom\^3\]:\s+([-\d.Ee+]+)")
RE_NATOMS = re.compile(r"- Atoms:\s+(\d+)")


def read(path) -> str:
    return Path(path).read_text(errors="replace")


def ended(text: str) -> bool:
    return "PROGRAM ENDED AT" in text


def scf_unconverged(text: str) -> bool:
    return ("SCF run NOT converged" in text) or ("outer SCF loop FAILED to converge" in text)


def opt_converged(text: str) -> bool:
    return "OPTIMIZATION COMPLETED" in text


def energies_ha(text: str) -> list[float]:
    return [float(x) for x in RE_ENERGY.findall(text)]


def final_energy_ev(text: str) -> float | None:
    e = energies_ha(text)
    return e[-1] * HARTREE_EV if e else None


def grid_counts(text: str) -> list[tuple[int, int, float]]:
    """(grid level, number of Gaussians, cutoff [Ry]) from the MULTIGRID INFO
    block (printed at PRINT_LEVEL MEDIUM). Only the first block is used."""
    rows, seen = [], set()
    for lvl, cnt, cut in RE_GRID.findall(text):
        lvl = int(lvl)
        if lvl in seen:
            break
        seen.add(lvl)
        rows.append((lvl, int(cnt), float(cut) * 2.0))  # a.u. (Ha) -> Ry
    return rows


def mpi_ranks(text: str) -> int | None:
    m = RE_NPROC.search(text)
    return int(m.group(1)) if m else None


def last_volume(text: str) -> float | None:
    v = RE_VOLUME.findall(text)
    return float(v[-1]) if v else None


def natoms(text: str) -> int | None:
    m = RE_NATOMS.search(text)
    return int(m.group(1)) if m else None


def find_output(run_dir: Path) -> Path | None:
    outs = sorted(run_dir.glob("*.inp.out")) or sorted(run_dir.glob("*.out"))
    outs = [o for o in outs if not o.name.startswith(("job", "slurm"))]
    return outs[0] if outs else None


def find_restart(run_dir: Path) -> Path | None:
    r = sorted(run_dir.glob("*-1.restart"))
    return r[0] if r else None
