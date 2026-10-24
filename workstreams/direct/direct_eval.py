"""Direct decision-observable measurement (EQS-style, Nakayama et al. arXiv 2408.03000) vs the
current best separate-overlap pipeline (W5). DEV ONLY.

M = sum_i w_i |psi(x_i)><psi(x_i)|, M = V diag(lam) V^dag. Circuit: U(x)|00>, then V^dag, measure; outcome j
scores lam_j. f = b + E[lam_J]. Checks: (1) exact equivalence; (2) noisy P(j) under Aer FakeFez on W5's subset
and pair; (3) adaptive Pocock rule on direct shots (T=200) noiseless and noisy, vs W5 OPTIMIZED.
"""
import sys
import numpy as np
from scipy.stats import norm

S = str(__import__('pathlib').Path(__file__).resolve().parents[1])
sys.path.insert(0, S + "/harness")
import qhack_harness as h
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeFez

d = h.load()
w, b, fmap = d["w"], d["b"], d["fmap"]
psi = lambda x: Statevector(fmap.assign_parameters(x)).data
M = sum(wi * np.outer(psi(xi), psi(xi).conj()) for wi, xi in zip(w, d["X_sv"]))
lam, V = np.linalg.eigh(M)
print("eigenvalues", lam.round(4), " trace(M) =", np.trace(M).real.round(12), " sum(w) =", w.sum().round(12))

# (1) exact equivalence on all dev points
P = np.array([np.abs(V.conj().T @ psi(x)) ** 2 for x in d["X_dev"]])          # (n_dev, 4) ideal outcome probs
f_dir = b + P @ lam
print(f"max |f_direct - f_exact| = {np.abs(f_dir - d['f_dev']).max():.1e}; sign disagreements "
      f"{(np.sign(f_dir) != np.sign(d['f_dev'])).sum()}")
var1 = P @ lam**2 - (P @ lam) ** 2                                          # per-shot variance of Z
K = d["K_dev"]
sep_neyman = (np.abs(w) * np.sqrt(K * (1 - K))).sum(1) ** 2                 # total-shot variance coefficient
print(f"ideal variance ratio (oracle-Neyman separate / direct): median {np.median(sep_neyman / var1):.2f}")

# (2) noisy outcome probabilities on W5's subset, W5's pair (4,5), opt level 3, seed 7
sub = np.load(S + "/w5_pipeline/probs.npz")["sub"]
be = FakeFez()
Vd = UnitaryGate(V.conj().T, label="Vdag")
def circ(x):
    qc = QuantumCircuit(2)
    qc.compose(fmap.assign_parameters(x), inplace=True)
    qc.append(Vd, [0, 1])
    qc.measure_all()
    return transpile(qc, be, initial_layout=[4, 5], optimization_level=3, seed_transpiler=7)
circs = [circ(d["X_dev"][i]) for i in sub]
print("direct circuit CZ counts:", np.bincount([c.count_ops().get("cz", 0) for c in circs]))
sim = AerSimulator.from_backend(be, seed_simulator=3)
SH = 200_000
res = sim.run(circs, shots=SH).result()
Pn = np.zeros((len(sub), 4))
for r, c in enumerate(circs):
    for key, v in res.get_counts(c).items():
        Pn[r, int(key.replace(" ", ""), 2)] += v / SH        # measure_all creg: bitstring 'q1q0' -> int = index
f_noisy = b + Pn @ lam
f_ex = d["f_dev"][sub]
near = d["near"][sub]
print(f"noisy direct: mean |f_noisy - f_exact| {np.abs(f_noisy - f_ex).mean():.4f}; "
      f"agreement with exact {np.mean(np.sign(f_noisy) == np.sign(f_ex)):.4f} (non-near "
      f"{np.mean((np.sign(f_noisy) == np.sign(f_ex))[~near]):.4f}); fitted slope "
      f"{np.polyfit(f_ex - b, f_noisy - b, 1).round(4)}")

# (3) adaptive Pocock on direct shots. Total shots per prediction on schedule 16*2^l, 9 looks (16..4096).
def adaptive_direct(Pmat, T=200, s0=16, looks=9, z=2.419, seed=0):
    rng = np.random.default_rng(seed)
    n_pts = len(Pmat)
    sched = s0 * 2 ** np.arange(looks)
    inc = np.diff(sched, prepend=0)
    cnt = np.zeros((T, n_pts, 4), np.int64)
    f_stop = np.zeros((T, n_pts)); s_stop = np.zeros((T, n_pts)); active = np.ones((T, n_pts), bool)
    for l, (n, add) in enumerate(zip(sched, inc)):
        cnt += np.stack([rng.multinomial(add, p, size=T) for p in Pmat], axis=1)
        pt = (cnt + 0.5) / (n + 2.0)                                       # pseudo-counts: never-zero variance
        mean_z = (cnt / n) @ lam
        var_z = pt @ lam**2 - (pt @ lam) ** 2
        f = b + mean_z
        stop = active & ((np.abs(f) >= z * np.sqrt(var_z / n)) | (l == looks - 1))
        f_stop[stop] = f[stop]; s_stop[stop] = n; active &= ~stop
    return f_stop, s_stop

def summarize(name, f_stop, s_stop, ref_sign, near):
    dis = (np.sign(f_stop) != ref_sign[None]).mean(0)
    print(f"{name:34s} agree {1 - dis.mean():.4f}  worst non-near {dis[~near].max():.3f}  "
          f"shots/pred mean {s_stop.mean():7.1f}  non-near mean {s_stop[:, ~near].mean():6.1f}  "
          f"p95 {np.percentile(s_stop, 95):6.0f}")

print("\nadaptive Pocock (z=2.419), T=200, W5 subset (155 pts: 150 far + 5 near-ties); reference = exact classifier")
ref = np.sign(f_ex)
summarize("direct, noiseless", *adaptive_direct(P[sub]), ref, near)
summarize("direct, FakeFez noise", *adaptive_direct(Pn), ref, near)
print("W5 OPTIMIZED (separate, noisy) for comparison: shots/pred mean 1859, non-near mean 698, worst non-near 0.010")
fs, ss = adaptive_direct(P)
dis = (np.sign(fs) != np.sign(d["f_dev"])[None]).mean(0)
print(f"\ndirect, noiseless, ALL 813 dev: agree {1 - dis.mean():.4f}, worst non-near {dis[~d['near']].max():.3f}, "
      f"shots/pred mean {ss.mean():.1f}, non-near mean {ss[:, ~d['near']].mean():.1f}; "
      f"dev accuracy {np.mean(np.where(fs > 0, d['classes'][1], d['classes'][0]) == d['y_dev'][None]):.4f}")
np.savez(S + "/direct/direct_probs.npz", lam=lam, V=V, P_dev=P, Pn=Pn, sub=sub)
