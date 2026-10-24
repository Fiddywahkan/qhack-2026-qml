import sys, numpy as np, json
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness')); sys.path.insert(0,".")
from qhack_harness import load
from circuits import *
from qiskit_ibm_runtime.fake_provider import FakeFez
from qiskit_aer import AerSimulator
be=FakeFez(); tg=be.target; d=load(2)
rng=np.random.default_rng(1)
samp=[(rng.integers(len(d["X_dev"])),rng.integers(len(d["X_sv"]))) for _ in range(3)]
def est_err(qc):
    ok=1.0; 
    for ins in qc.data:
        n=ins.operation.name
        if n in("barrier",): continue
        q=tuple(qc.find_bit(b).index for b in ins.qubits)
        if n=="measure": e=tg["measure"][q].error
        elif n in tg.operation_names and q in tg[n]: e=tg[n][q].error
        else: continue
        ok*=1-(e or 0)
    gate=1-ok
    dur=qc.estimate_duration(tg); qs=sorted({qc.find_bit(b).index for ins in qc.data for b in ins.qubits if ins.operation.name!="barrier"})
    dec=0
    for q in qs:
        p=be.qubit_properties(q); dec+=dur*(0.5/p.t1+0.5/p.t2)
    return gate, dec
rng=np.random.default_rng(2)

import pickle
res=pickle.load(open('layout.pkl','rb')); best=res[0][3]; medp=res[len(res)//2][3]
sim=AerSimulator.from_backend(be)
idx=[(rng.integers(len(d["X_dev"])),rng.integers(len(d["X_sv"]))) for _ in range(20)]
S=4096
for pname,pair in (("best",best),("median",medp)):
  for var,lvl in (("original",0),("original",1),("original",3),("merged",0),("merged",3)):
    errs=[];ks=[]
    for n,(i,j) in enumerate(idx):
        x,xi=d["X_dev"][i],d["X_sv"][j]
        kex=p00(merged(x,xi,2))
        qc=build_transpiled(x,xi,be,pair,lvl,var,d["fmap"],2)
        c=sim.run(qc,shots=S,seed_simulator=100+n).result().get_counts()
        kn=c.get("00",0)/S; errs.append(kn-kex)
    errs=np.array(errs)
    print(pname,pair,var,"lvl",lvl,"mean|dk| %.4f  mean bias %.4f"%(np.abs(errs).mean(),errs.mean()))
print("shot-noise floor approx mean|dk| ~ %.4f"%(np.mean([np.sqrt(0.25/S)*0.8])))
