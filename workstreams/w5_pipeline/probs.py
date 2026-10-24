"""Compute noisy P(00) (FakeFez, Aer density_matrix, 2e5 shots -> sd~1e-3) for every (dev subset x SV) circuit,
for (a) W3 circuit pair (4,5) and (b) W4 parametrized template layout [0,1]. Plus calibration/identity circuits."""
import sys, time, numpy as np
S=str(__import__('pathlib').Path(__file__).resolve().parents[1])
sys.path[:0]=[S+'/harness', S+'/w3_circuit']
from multiprocessing import Pool
from qhack_harness import load
from recommended import kernel_circuit
from circuits import merged, p00
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeFez
SH=200_000
d=load(); be=FakeFez(); sim=AerSimulator.from_backend(be, method='density_matrix')
rng=np.random.default_rng(0)   # subset seed 0
near=np.where(d['near'])[0]; far=np.setdiff1d(np.arange(len(d['X_dev'])),near)
sub=np.sort(np.concatenate([rng.choice(far,150,replace=False),near]))
# template
x=ParameterVector('x',2); y=ParameterVector('y',2); qc=QuantumCircuit(2)
qc.compose(d['fmap'].assign_parameters(list(x)),inplace=True); qc.compose(d['fmap'].assign_parameters(list(y)).inverse(),inplace=True); qc.measure_all()
T=transpile(qc,be,initial_layout=[0,1],optimization_level=3,seed_transpiler=11)
def tb(a,b): return T.assign_parameters({p:(a if p.name.startswith('x[') else b)[int(p.name[2:-1])] for p in T.parameters})
def w3(a,b): return kernel_circuit(a,b,be,pair=(4,5))
def job(chunk):
    cs=[(tb(a,b) if kind=='t' else w3(a,b)) for kind,a,b in chunk]
    r=sim.run(cs,shots=SH,seed_simulator=5,max_parallel_experiments=16).result()
    return [(r.get_counts(i).get('00',0)/SH, c.count_ops().get('cz',0), c.depth()) for i,c in enumerate(cs)]
if __name__=='__main__':
    t0=time.time()
    # calibration pairs from TRAIN (non-SV train point vs SV): choose 8 with exact k farthest from 1/4 (highest |k-1/4|)
    Xt=d['X_train']; nonsv=np.setdiff1d(np.arange(len(Xt)),d['sv'])
    cand=[(i,j) for i in nonsv for j in range(20)]
    kc=np.array([p00(merged(Xt[i],d['X_sv'][j])) for i,j in cand])
    order=np.argsort(-np.abs(kc-0.25))[:8]; cal=[cand[o] for o in order]; kcal=kc[order]
    jobs=[]
    for kind in 'tw':
        for i in sub:
            for j in range(20): jobs.append((kind,d['X_dev'][i],d['X_sv'][j]))
        for j in range(4): jobs.append((kind,d['X_sv'][j],d['X_sv'][j]))          # identity
        for (i,j) in cal: jobs.append((kind,Xt[i],d['X_sv'][j]))                  # cal pairs
    chunks=[jobs[i:i+400] for i in range(0,len(jobs),400)]; print('njobs',len(jobs),flush=True)
    out=[job(ch) for ch in chunks]
    out=np.array([o for ch in out for o in ch]); n=len(sub)*20+4+8
    res={}
    for q,kind in enumerate('tw'):
        o=out[q*n:(q+1)*n]
        res[kind+'_P']=o[:len(sub)*20,0].reshape(len(sub),20); res[kind+'_id']=o[len(sub)*20:len(sub)*20+4,0]
        res[kind+'_cal']=o[len(sub)*20+4:,0]; res[kind+'_cz']=o[:,1]; res[kind+'_depth']=o[:,2]
        print(kind,'mean cz dev',o[:len(sub)*20,1].mean(),'identity cz/depth',o[len(sub)*20:len(sub)*20+4,1:].tolist())
    np.savez('probs.npz',sub=sub,kcal=kcal,**res)
    print('time',time.time()-t0)
