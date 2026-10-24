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
edges=sorted({tuple(sorted(e)) for e in be.coupling_map.get_edges()})
res=[]
for e in edges:
    g=[];dc=[]
    for i,j in samp:
        qc=build_transpiled(d["X_dev"][i],d["X_sv"][j],be,e,3,"merged",None,2)
        a,b=est_err(qc); g.append(a);dc.append(b)
    res.append((np.mean(g)+np.mean(dc),np.mean(g),np.mean(dc),e))
res.sort()
import pickle;pickle.dump(res,open('layout.pkl','wb'))
print("n edges",len(res))
print("best5",[(round(r[0],4),r[3]) for r in res[:5]])
med=res[len(res)//2]; print("median",round(med[0],4),med[3],"gate",round(med[1],4),"dec",round(med[2],4))
print("worst",round(res[-1][0],4),res[-1][3])
print("best gate/dec",round(res[0][1],4),round(res[0][2],4))
best=res[0][3]; medp=med[3]
print("BEST",best,"MED",medp)
