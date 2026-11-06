import re,glob,os,sys
B=sys.argv[1]
H=27211.386245988
def seg(f):
    t=open(f,errors='ignore').read()
    i=t.rfind('PROGRAM STARTED AT')
    return t[i:]
def nat(s):
    m=re.search(r'- Atoms:\s+(\d+)',s); return int(m.group(1)) if m else None
def energy(s):
    m=re.findall(r'ENERGY\| Total FORCE_EVAL.*?(-?\d+\.\d+)\s*$',s,re.M); return float(m[-1]) if m else None
def vol(s):
    m=re.findall(r'CELL\| Volume \[angstrom\^3\]:\s+(\d+\.\d+)',s); return float(m[-1]) if m else None
def grids(s):
    return [int(x) for x in re.findall(r'count for grid\s+\d+:\s+(\d+)',s)]
for fn in ['PBE','PBEsol','revPBE','BLYP']:
    print('\n=====',fn)
    for kind,pre in (('cutoff','cutoff'),('relcutoff','relcutoff')):
        rows=[]
        for d in glob.glob(f'{B}/{fn}/grid/{pre}*'):
            n=int(re.sub(r'\D','',os.path.basename(d)))
            if kind=='cutoff' and 'rel' in os.path.basename(d): continue
            f=glob.glob(d+'/*.out')
            if not f: continue
            s=seg(f[0]); rows.append((n,energy(s),grids(s),nat(s)))
        rows.sort()
        if not rows: continue
        N=rows[0][3] or 1; ref=rows[-1][1]
        print(f'-- {kind} (N={N}; ref=last)  n | dE[meV/at] | Gaussians per grid')
        for n,e,g,_ in rows:
            print(f'  {n:5d} {(e-ref)*H/N:9.3f}  {g}')
    print('-- kmesh CELL_OPT (converged only)')
    rows=[]
    for d in glob.glob(f'{B}/{fn}/kmesh/k*'):
        if 'failed' in d: continue
        f=glob.glob(d+'/*.out')
        if not f: continue
        s=seg(f[0]); ok='OPTIMIZATION COMPLETED' in s
        rows.append((os.path.basename(d),ok,energy(s),vol(s),nat(s)))
    ref=[r for r in rows if r[1]]
    if not ref: continue
    ref=sorted(ref,key=lambda r:int(r[0][1:].split('_')[0]))[-1]
    for k,ok,e,v,N in sorted(rows,key=lambda r:int(r[0][1:].split('_')[0])):
        if ok: print(f'  {k:8s} OK   dE={(e-ref[2])*H/N:8.3f} meV/at  V={v:.3f} dV={100*(v-ref[3])/ref[3]:6.2f}%')
        else: print(f'  {k:8s} --   (não fechou)')
    print('  ref =',ref[0])
