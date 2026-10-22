#!/usr/bin/env python3
"""Plot xz and yz projections of one or more xyz files side by side (quick visual check).

Usage: preview_xyz.py out.png a.xyz [b.xyz ...]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = {"Au": "#d4af37", "C": "#555555", "N": "#2255cc", "H": "#bbbbbb"}
SIZE = {"Au": 220, "C": 90, "N": 100, "H": 30}


def read(path):
    L = open(path).read().split("\n")
    n = int(L[0])
    el = [l.split()[0] for l in L[2:2 + n]]
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return el, xyz


files = sys.argv[2:]
fig, axs = plt.subplots(2, len(files), figsize=(3.2 * len(files), 8), squeeze=False)
for j, f in enumerate(files):
    el, X = read(f)
    for i, (a, b, lab) in enumerate(((0, 2, "x"), (1, 2, "y"))):
        ax = axs[i, j]
        order = np.argsort(X[:, 1 - i])
        for k in order:
            ax.scatter(X[k, a], X[k, b], s=SIZE[el[k]], c=COL[el[k]], edgecolors="k", linewidths=0.3)
        ax.set_aspect("equal"); ax.set_xlabel(f"{lab} (A)"); ax.set_ylabel("z (A)")
        if i == 0:
            ax.set_title(f.split("/")[-1].replace(".xyz", ""), fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=110)
