"""analyze_sigma_etot.py [ROOT] [SUBDIR] -- E_free, -TS, E0=E-corr/2 vs sigma per functional (last PROGRAM STARTED segment).
ROOT default '.' (benchmark dir); SUBDIR default 'etot' (use 'etot_relaxed' with layout F/sigmaX -> pass SUBDIR=. and ROOT=etot_relaxed)."""
import re,sys,os
ROOT=sys.argv[1] if len(sys.argv)>1 else "."
SUB=sys.argv[2] if len(sys.argv)>2 else "etot"
N=24; H=27211.386
S=[0.01,0.1,0.5,1,2,3,4,5,6]
for f in ["PBE","PBEsol","revPBE","BLYP"]:
    print(f"\n== {f}")
    rows=[]
    for s in S:
        n=f"sigma{s:g}"; o=os.path.join(ROOT,f,SUB,n,f"Ag2SeO3_{f}_{n}.out") if SUB!="." else os.path.join(ROOT,f,n,f"Ag2SeO3_{f}_{n}.out")
        if not os.path.exists(o): print(n,"sem .out"); continue
        t=open(o).read(); t=t[t.rfind("PROGRAM STARTED AT"):]
        ok="OPTIMIZATION COMPLETED" in t; nc="SCF run NOT converged" in t
        E=float(re.findall(r"ENERGY\| Total FORCE_EVAL.*?\]\s+(-?[\d.]+)",t)[-1])
        c=float(re.findall(r"Smearing free energy correction:\s+(-?[\d.]+)",t)[-1])
        mpi=re.search(r"message passing processes\s+(\d+)",t)
        stp=len(re.findall(r"OPTIMIZATION STEP",t))
        g=re.findall(r"Max. gradient\s+=\s+(\S+)",t)
        rows.append((s,E,c,E-c/2,ok,nc,stp,g[-1] if g else "?",mpi.group(1) if mpi else "?"))
    if not rows: continue
    r0=rows[0][3]
    print("sigma  E_free(Ha)   -TS(meV/at)  E0-E0ref(meV/at)  opt scf_ok steps maxgrad mpi")
    for s,E,c,E0,ok,nc,stp,g,m in rows:
        print(f"{s:<6g} {E:.6f} {c*H/N:9.1f} {(E0-r0)*H/N:10.2f}   {ok} {not nc} {stp} {g} {m}")
