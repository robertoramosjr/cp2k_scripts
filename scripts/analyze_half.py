#!/usr/bin/env python3
"""Read a half-bridge GEO_OPT trajectory (<project>-pos-1.xyz): per frame energy, S/N-bottom-Au
and closest molecule-top-tip distances, split by element (Au 1-14 = bottom tip, 15-28 = top tip).

usage: analyze_half.py <project>-pos-1.xyz [step]
"""
import sys

import numpy as np

path, step = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10
L = open(path).read().split("\n")
n = int(L[0]); b = n + 2
frames = [L[i:i + b] for i in range(0, len(L) - b + 1, b) if L[i].strip()]
print(f"{len(frames)} frames, {n} atoms; columns: frame  E[Eh]  Au_bot-anchor  top-tip nearest (C | H | N,S)")
for k in sorted(set(range(0, len(frames), step)) | {len(frames) - 1}):
    F = frames[k]
    E = float(F[1].split("E =")[1])
    el = [l.split()[0] for l in F[2:2 + n]]
    P = np.array([[*map(float, l.split()[1:4])] for l in F[2:2 + n]])
    bot, top = P[:14], P[14:28]
    mi = [i for i in range(28, n)]
    anc = [i for i in mi if el[i] in ("S",)] or [i for i in mi if el[i] == "N" and sum(np.linalg.norm(P[j] - P[i]) < 1.2 for j in mi if el[j] == "H") == 2]
    da = min(np.linalg.norm(bot - P[i], axis=1).min() for i in anc)
    nd = lambda syms: min((np.linalg.norm(top - P[i], axis=1).min() for i in mi if el[i] in syms), default=float("nan"))
    print(f"{k + 1:5d} {E:16.6f}  {da:5.2f}   {nd('C'):5.2f} | {nd('H'):5.2f} | {nd('NS'):5.2f}")
