#!/usr/bin/env python3
"""
parse_cutoff_convergence.py

Parses CP2K output files produced by cutoff_convergence.py, computes the
total-energy change (meV/atom) relative to the highest CUTOFF tested, plots
E vs CUTOFF, and suggests a converged (cost-efficient) cutoff -- the CP2K
equivalent of a VASP ENCUT convergence plot.

Usage:
    python ~/scripts/parse_cutoff_convergence.py \\
        --scan-dir cutoff_convergence --n-atoms 6 --threshold-mev-atom 1.0
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ENERGY_PATTERN = re.compile(r"(?m)^\s*ENERGY\|\s+Total FORCE_EVAL.*\s(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*$")
# Anchored at line start so it does NOT match "REL_CUTOFF   60"
CUTOFF_PATTERN = re.compile(r"(?m)^\s*CUTOFF\s+(\d+)")

HARTREE_TO_MEV = 27211.386245988


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse and plot a CP2K MGRID CUTOFF convergence scan."
    )
    parser.add_argument("--scan-dir", type=str, default="cutoff_convergence",
                         help="Root directory produced by cutoff_convergence.py")
    parser.add_argument("--n-atoms", type=int, required=True,
                         help="Number of atoms in the cell (for the meV/atom criterion).")
    parser.add_argument("--threshold-mev-atom", type=float, default=1.0,
                         help="Convergence threshold in meV/atom relative to the highest cutoff.")
    parser.add_argument("--output-plot", type=str, default="cutoff_convergence.png")
    return parser.parse_args()


def extract_energy_hartree(out_file: Path):
    text = out_file.read_text(errors="ignore")
    matches = ENERGY_PATTERN.findall(text)
    return float(matches[-1]) if matches else None


def scf_converged(out_file: Path):
    # CP2K does NOT abort on SCF non-convergence -- it still prints an
    # ENERGY| line with whatever unconverged density it has. Require the
    # explicit outer-SCF convergence marker, or a bad point would silently
    # end up in the convergence curve with no indication anything was wrong.
    text = out_file.read_text(errors="ignore")
    return "outer SCF loop converged" in text


def main():
    args = parse_args()
    scan_root = Path(args.scan_dir)
    if not scan_root.exists():
        raise SystemExit(f"ERROR: scan directory not found: {scan_root}")

    results = []  # list of (cutoff, energy_hartree)
    for run_dir in sorted(scan_root.iterdir()):
        if not run_dir.is_dir():
            continue

        out_files = list(run_dir.glob("*.out")) + list(run_dir.glob("*.inp.out"))
        if not out_files:
            print(f"WARNING: no CP2K output found in {run_dir}, skipping.")
            continue

        if not scf_converged(out_files[0]):
            print(f"WARNING: SCF did NOT converge in {out_files[0]}, skipping.")
            continue

        energy_ha = extract_energy_hartree(out_files[0])
        if energy_ha is None:
            print(f"WARNING: could not parse total energy in {out_files[0]}, skipping.")
            continue

        inp_files = list(run_dir.glob("*.inp"))
        cutoff = None
        if inp_files:
            m = CUTOFF_PATTERN.search(inp_files[0].read_text())
            if m:
                cutoff = int(m.group(1))
        if cutoff is None:
            print(f"WARNING: could not read CUTOFF from input in {run_dir}, skipping.")
            continue

        results.append((cutoff, energy_ha))

    if not results:
        raise SystemExit("No results parsed. Check --scan-dir and that CP2K runs completed.")

    results.sort(key=lambda r: r[0])

    cutoffs = [r[0] for r in results]
    energies_ha = [r[1] for r in results]

    reference = energies_ha[-1]  # highest cutoff tested = reference
    delta_mev_atom = [(e - reference) * HARTREE_TO_MEV / args.n_atoms for e in energies_ha]

    print(f"{'CUTOFF (Ry)':<14}{'E (Ha)':>18}{'dE (meV/atom)':>18}")
    for cutoff, e_ha, d in zip(cutoffs, energies_ha, delta_mev_atom):
        print(f"{cutoff:<14}{e_ha:>18.8f}{d:>18.4f}")

    recommended = cutoffs[-1]
    for cutoff, d in zip(cutoffs, delta_mev_atom):
        if abs(d) <= args.threshold_mev_atom:
            recommended = cutoff
            break

    print(f"\nRecommended CUTOFF (first one within {args.threshold_mev_atom} meV/atom "
          f"of the highest cutoff tested): {recommended} Ry")
    print("Use this value as --cutoff in kpoint_convergence.py for the next step.")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(cutoffs, delta_mev_atom, marker="o")
    ax.axhline(args.threshold_mev_atom, color="gray", linestyle="--", linewidth=1,
               label=f"±{args.threshold_mev_atom} meV/atom threshold")
    ax.axhline(-args.threshold_mev_atom, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("CUTOFF (Ry)")
    ax.set_ylabel("ΔE relative to highest cutoff (meV/atom)")
    ax.set_title("MGRID cutoff convergence")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_plot, dpi=150)
    print(f"\nPlot saved to: {args.output_plot}")


if __name__ == "__main__":
    main()