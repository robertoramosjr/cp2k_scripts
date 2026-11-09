#!/usr/bin/env python3
"""check_cellopt.py -- verify finished CELL_OPT runs: OPTIMIZATION COMPLETED, no SCF failure, and
|a|,|b|,|c| of the final cell vs CELL_REF (fixed-grid reference). Last 'PROGRAM STARTED AT' segment only.
Usage: check_cellopt.py <inp>.out ..."""
import re, sys
for out in sys.argv[1:]:
    t = open(out).read(); t = t[t.rfind("PROGRAM STARTED AT"):]
    ok = "OPTIMIZATION COMPLETED" in t and "SCF run NOT converged" not in t
    inp = re.sub(r"\.out$", ".inp", out)
    ref = re.search(r"&CELL_REF\n(.*?)&END CELL_REF", open(inp).read(), re.S).group(1)
    nref = [sum(float(x)**2 for x in re.search(rf"^\s*{k} \[angstrom\]\s+(.*)$", ref, re.M).group(1).split())**.5 for k in "ABC"]
    fin = [float(x) for x in re.findall(r"\|[abc]\| =\s+([\d.]+)", t)[-3:]]
    flag = "OK" if all(f <= r for f, r in zip(fin, nref)) else "EXCEEDS CELL_REF -> redo with make_cellopt_from_restart.py"
    print(f"{out}\n  completed={ok}  final |a|,|b|,|c|={fin}  CELL_REF={[round(x,3) for x in nref]}  {flag}")
