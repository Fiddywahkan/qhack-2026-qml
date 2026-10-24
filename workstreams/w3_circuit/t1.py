import sys, json, numpy as np
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness'))
sys.path.insert(0,".")
from qhack_harness import load
from circuits import *
from qiskit_ibm_runtime.fake_provider import FakeFez
be=FakeFez(); d=load(2)
rng=np.random.default_rng(0)
# equivalence check
err=0
for r in (1,2,3):
    dd=load(r) if r!=2 else d
    for _ in range(30):
        x=dd["X_dev"][rng.integers(len(dd["X_dev"]))]; xi=dd["X_sv"][rng.integers(len(dd["X_sv"]))]
        err=max(err,abs(p00(original(dd["fmap"],x,xi))-p00(merged(x,xi,r))))
print("max |k_orig-k_merged| over reps1-3, 90 pairs:",err)
pair=(0,1)
for r in (2,):
  idx=[(rng.integers(len(d["X_dev"])),rng.integers(len(d["X_sv"]))) for _ in range(20)]
  for var in ("original","merged"):
    for lvl in range(4):
        st=[]
        for i,j in idx:
            qc=build_transpiled(d["X_dev"][i],d["X_sv"][j],be,pair,lvl,var,d["fmap"],r)
            ops=qc.count_ops(); cz=ops.get("cz",0)
            n1=sum(v for k,v in ops.items() if k in("rz","sx","x"))
            dur=qc.estimate_duration(be.target)*1e0 if hasattr(qc,"estimate_duration") else np.nan
            st.append((cz,n1,qc.depth(),dur))
        a=np.array(st); print(var,lvl,"CZ %.2f 1q %.1f depth %.1f dur_us %.2f"%(a[:,0].mean(),a[:,1].mean(),a[:,2].mean(),a[:,3].mean()*1e6))
