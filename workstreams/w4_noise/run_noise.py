import sys, json, time, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness'))
sys.path.insert(0, '.')
from qhack_harness import load
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeFez
SHOTS, NSUB, SEED = 20000, 150, 7
d = load(); fmap = d['fmap']
be = FakeFez(); sim = AerSimulator.from_backend(be)
edge = [e for e in be.coupling_map.get_edges() if e[0] < e[1]][0]
layout = list(edge); print('layout', layout)
x = ParameterVector('x', 2); y = ParameterVector('y', 2)
qc = QuantumCircuit(2)
qc.compose(fmap.assign_parameters(list(x)), inplace=True)
qc.compose(fmap.assign_parameters(list(y)).inverse(), inplace=True)
qc.measure_all()
t = transpile(qc, be, initial_layout=layout, optimization_level=3, seed_transpiler=11)
print('transpiled ops', dict(t.count_ops()), 'depth', t.depth())
# calibration circuits (readout)
cal = []
for p in ['00','01','10','11']:
    c = QuantumCircuit(2)
    if p[1]=='1': c.x(0)
    if p[0]=='1': c.x(1)
    c.measure_all(); cal.append(transpile(c, be, initial_layout=layout, optimization_level=0, seed_transpiler=11))
labels=['00','01','10','11']
def dist(counts): return np.array([counts.get(l,0) for l in labels])/SHOTS
def run(circs, seed):
    r = sim.run(circs, shots=SHOTS, seed_simulator=seed).result()
    return np.array([dist(r.get_counts(i)) for i in range(len(circs))])
def bind(xa, ya):
    ps = list(t.parameters); m={}
    for p in ps:
        m[p] = (xa if p.name.startswith('x') else ya)[int(p.name[2:-1])] if False else None
    return m
def bound(xa, ya):
    return t.assign_parameters({p: (xa if p.name.startswith('x[') else ya)[int(p.name[2:-1])] for p in t.parameters})
rng = np.random.default_rng(SEED)
nearidx = np.where(d['near'])[0]
far = np.setdiff1d(np.arange(len(d['X_dev'])), nearidx)
sub = np.sort(np.concatenate([rng.choice(far, NSUB, replace=False), nearidx]))
print('n dev subset', len(sub), 'near', len(nearidx))
pairs = [(i,j) for i in sub for j in range(20)]
t0=time.time()
Pcal = run(cal, 1)
A = np.array(Pcal).T  # A[meas,prep]
print('A diag', np.diag(A))
Pd = run([bound(d['X_dev'][i], d['X_sv'][j]) for i,j in pairs], 2).reshape(len(sub),20,4)
print('dev run s', time.time()-t0)
# calibration set for affine fit: training points vs SVs (training data, not dev): 40 random pairs + 4 identity circuits
Xt = d['X_train']; ci = rng.choice(len(Xt), 40, replace=False)
cal_pairs = [(Xt[a], d['X_sv'][j]) for a,j in zip(ci, rng.integers(0,20,40))]
Pc = run([bound(a,b_) for a,b_ in cal_pairs], 3)
kc_exact = np.array([d['qkernel'].evaluate(a[None], b_[None])[0,0] for a,b_ in cal_pairs])
idc = [(d['X_sv'][j], d['X_sv'][j]) for j in range(4)]
Pid = run([bound(a,b_) for a,b_ in idc], 4)
np.savez('raw.npz', Pd=Pd, sub=sub, A=A, Pc=Pc, kc_exact=kc_exact, Pid=Pid)
