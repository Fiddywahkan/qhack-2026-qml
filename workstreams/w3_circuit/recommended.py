"""Recommended hardware kernel circuit: fused ('merged') ZZ map, FakeFez/ibm_fez, level 3, seed 7.
Usage: qc = kernel_circuit(x, xi, backend, pair=(4,5)); (pair chosen on FakeFez calibration: re-pick on live ibm_fez calibration!)"""
from circuits import build_transpiled, merged, p00
def kernel_circuit(x, xi, backend, pair=(4, 5), level=3, seed=7, reps=2):
    return build_transpiled(x, xi, backend, pair, level, "merged", None, reps, seed)
