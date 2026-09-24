#!/usr/bin/env python3
"""
extract_relaxed_structure.py

Extracts the final relaxed &CELL and &COORD from a CP2K .restart file
(written after CELL_OPT/GEO_OPT) and writes it as a CIF via pymatgen --
ready to feed back into input_generator.py (--cif ...) for the next
pipeline step (e.g. taking the relaxed structure from the winning k-mesh
in kmesh_stress_convergence.py into the production Etot/GEO_OPT run).

Usage:
    python ~/scripts/extract_relaxed_structure.py \\
        --restart kmesh_stress_convergence/k6-6-6/stress_k6-6-6-1.restart \\
        --output TiO2_anatase_relaxed.cif
"""

import argparse
import re
import sys
from pathlib import Path

try:
    from pymatgen.core import Structure, Lattice
except ImportError:
    sys.exit("ERROR: pymatgen is required. Install with: python3 -m pip install pymatgen")

CELL_VECTOR_PATTERN = re.compile(
    r"^\s*([ABC])\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)", re.MULTILINE
)
COORD_BLOCK_PATTERN = re.compile(r"&COORD(.*?)&END COORD", re.DOTALL)
COORD_LINE_PATTERN = re.compile(
    r"^\s*([A-Za-z]{1,2})\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)", re.MULTILINE
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a CP2K .restart file's final geometry into a CIF."
    )
    parser.add_argument("--restart", type=str, required=True,
                         help="Path to the CP2K *-1.restart file.")
    parser.add_argument("--output", type=str, default=None,
                         help="Output CIF path. Defaults to <restart_stem>.cif")
    return parser.parse_args()


def parse_cell(text: str):
    vectors = {}
    for m in CELL_VECTOR_PATTERN.finditer(text):
        label, x, y, z = m.groups()
        vectors[label] = (float(x), float(y), float(z))  # last occurrence wins
    if not {"A", "B", "C"} <= vectors.keys():
        sys.exit("ERROR: could not find a complete &CELL (A/B/C) block in the restart file.")
    return Lattice([vectors["A"], vectors["B"], vectors["C"]])


def parse_coords(text: str):
    coord_matches = COORD_BLOCK_PATTERN.findall(text)
    if not coord_matches:
        sys.exit("ERROR: could not find a &COORD block in the restart file.")
    coord_block = coord_matches[-1]  # last occurrence = final geometry
    is_scaled = "SCALED" in coord_block.upper()

    species, coords = [], []
    for m in COORD_LINE_PATTERN.finditer(coord_block):
        symbol, x, y, z = m.groups()
        species.append(symbol)
        coords.append((float(x), float(y), float(z)))

    if not species:
        sys.exit("ERROR: &COORD block found but no atom lines could be parsed.")
    return species, coords, is_scaled


def main():
    args = parse_args()
    restart_path = Path(args.restart)
    if not restart_path.exists():
        sys.exit(f"ERROR: restart file not found: {restart_path}")

    text = restart_path.read_text(errors="ignore")
    lattice = parse_cell(text)
    species, coords, is_scaled = parse_coords(text)

    structure = Structure(
        lattice, species, coords, coords_are_cartesian=not is_scaled
    )

    output_path = Path(args.output) if args.output else restart_path.with_suffix(".cif")
    structure.to(filename=str(output_path))

    print(f"Parsed {len(species)} atoms, coords_type={'SCALED' if is_scaled else 'CARTESIAN'}")
    print(f"Lattice: a={lattice.a:.5f} b={lattice.b:.5f} c={lattice.c:.5f} "
          f"(alpha={lattice.alpha:.2f} beta={lattice.beta:.2f} gamma={lattice.gamma:.2f})")
    print(f"CIF written to: {output_path}")
    print(f"\nFeed it into the next step with, e.g.:")
    print(f"  python ~/scripts/input_generator.py --cif {output_path} "
          f"--project-name TiO2_etot --run-type GEO_OPT --cutoff <CUTOFF> --kpoints <KX> <KY> <KZ>")


if __name__ == "__main__":
    main()