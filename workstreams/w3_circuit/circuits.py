"""Kernel circuit builders (original and analytically merged) + transpile helper. Seeds explicit."""
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import zz_feature_map
from qiskit.quantum_info import Statevector

def original(fmap, x, xi):
    qc = fmap.assign_parameters(x).compose(fmap.assign_parameters(xi).inverse())
    qc.measure_all()
    return qc

def _diag(qc, x):
    # one ZZ-map diagonal layer: P(2x_j) each, then CX P(2(pi-x0)(pi-x1)) CX
    return 2*x[0], 2*x[1], 2*(np.pi-x[0])*(np.pi-x[1])

def merged(x, xi, reps=2):
    """U(xi)^dag U(x) with the adjacent diagonal layers D(xi)^dag D(x) fused into one
    (P(2(x_j-xi_j)) and one RZZ-type phase). Equivalent to zz_feature_map(2,reps) kernel up to global phase.
    CX-P(t)-CX == RZZ(t) up to global phase."""
    qc = QuantumCircuit(2)
    def Ulayer_front(v):   # H then D(v)
        qc.h([0,1]); a,b,t = _diag(qc, v); qc.p(a,0); qc.p(b,1); qc.rzz(t,0,1)
    def Ulayer_front_inv(v):  # (H D(v))^dag = D(v)^dag H
        a,b,t = _diag(qc, v); qc.rzz(-t,0,1); qc.p(-a,0); qc.p(-b,1); qc.h([0,1])
    # U(x) = D_R H ... D_1 H  (first r-1 reps full, last rep: H only, D fused)
    for _ in range(reps-1): Ulayer_front(x)
    qc.h([0,1])
    # fused middle diagonal D(x) D(xi)^dag
    a,b,t = _diag(qc,x); a2,b2,t2 = _diag(qc,xi)
    qc.p(a-a2,0); qc.p(b-b2,1); qc.rzz(t-t2,0,1)
    qc.h([0,1])
    for _ in range(reps-1): Ulayer_front_inv(xi)
    qc.measure_all()
    return qc

def p00(qc):
    q = qc.remove_final_measurements(inplace=False)
    return float(abs(Statevector(q).data[0])**2)

def build_transpiled(x, xi, backend, pair, level=3, variant="merged", fmap=None, reps=2, seed=7):
    qc = merged(x, xi, reps) if variant == "merged" else original(fmap, x, xi)
    return transpile(qc, backend, initial_layout=list(pair), optimization_level=level, seed_transpiler=seed)
