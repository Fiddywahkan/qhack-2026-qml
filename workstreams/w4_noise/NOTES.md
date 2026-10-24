# W4 noise robustness (FakeFez, Aer, dev only)
Setup: FakeFez noise model, layout [0,1], transpile opt-level 3, seed_transpiler=11 (one parametrized circuit: 30 rz, 16 sx, 6 cz, depth 40), 20000 shots/circuit (shot sd ~0.0035), seeds 1-4 sim, 7 subset. Subset = 150 random non-near + 5 near dev points (155) x 20 SVs = 3100 circuits.

| setting | extra circuits (shots) | mean abs(k-k_exact) | max | agree w/ exact (155 pts) | dev acc |
|---|---|---|---|---|---|
| unmitigated | 0 | 0.0204 | 0.086 | 0.994 | 0.987 |
| readout (RO) mitigation | 4 (4x20k) | 0.0153 | 0.067 | 0.994 | 0.987 |
| depol rescale, identity circuits only | 4 | 0.0033 | 0.016 | 0.994 | 0.987 |
| RO + depol rescale (identity) | 4+4 | 0.0028 | 0.013 | 0.994 | 0.987 |
| affine fit, 40 cal pairs (no RO) | 40 | 0.0033 | 0.016 | 0.994 | 0.987 |
| RO + affine fit, 40 pairs | 4+40 | 0.0027 | 0.015 | 0.994 | 0.987 |

Exact classifier dev acc on the subset: 0.981. Shot noise floor in k is ~0.0035/circuit, so residuals after depol rescale (~0.003) are at the noise floor.
The single disagreement (all rows) is a near-boundary point with |f_exact| = 0.0009; no non-near point flips (far agreement 1.000).

## Findings
- Bias is almost exactly affine: k_noisy = 0.888 k + 0.028 (residual sd 0.0037, corr 0.9998; quadratic term -7e-4). After RO mitigation: 0.916 k + 0.020. Identity circuits give k(x,x)=0.916 (RO-mitigated 0.937), consistent with a = (k_id - 1/4)/(3/4) = 0.888.
- KEY: SVC dual constraint gives sum(w) = 7.8e-16 = 0, so the offset c cancels exactly in f. f_noisy = a(f - b) + b, i.e. the sign threshold shifts by b(1-a) = -0.0141 (b = -0.1255). Only points with |f| below ~0.014 can flip. Full dev (813 pts, model-based from a=0.888, script extra.py): 2 flips, both in the 5 near-tie points. So on this task unmitigated noise barely changes the classifier; the "needs mitigation" effect is limited to near-ties.
- Correct mitigated rule: f_corr = sum w_i (k_i - c)/a + b. Because sum(w)=0 this equals (f_noisy - b)/a + b: only a is needed (c irrelevant). a from identity circuits assuming c=(1-a)/4 matched the 40-pair fit to 0.001.
- Mitigation gains little for sign agreement here; main value is removing the threshold shift for near-ties. Cost of 1/a rescale: shot-variance inflated by 1/a^2 ~ 1.2-1.3x (i.e. ~20% more shots) -- not free in the adaptive regime.
- Robustness: RO+fit with 5, 10, 20 random cal pairs gives the same agreement (0.994) in all 200 resamples.

## Recommendation for hardware
Cheapest: 4 identity circuits k(x_sv,x_sv)=1 (use distinct SVs, ~2-4k shots each, <1 QPU-second at Heron rates) -> a = (mean k_id - 1/4)/(3/4); apply f = (f_noisy - b)/a + b (c ignored since sum w = 0; verify sum w ~ 0 on the actual model). Optionally add 4 RO calibration circuits (prep 00/01/10/11) and apply assignment-matrix inverse first (mitigation.py) for a slightly better model; or enable runtime TREX/measurement twirling (resilience/twirling options, e.g. twirling.enable_measure, resilience.measure_mitigation; from memory of qiskit-ibm-runtime docs, not run) which yields the same effect at no extra circuits of mine. Identity circuits must be transpiled identically (same layout, depth; note an x,x circuit may be compiled shallower by opt-level 3 -> build it with the same parametrized template and assign equal params, and check transpiled depth). Dynamical decoupling: circuits are ~40-deep, qubit idle time is minimal, expected irrelevant. Avoid ZNE/PEC: 3x+ circuit cost for no needed gain.
## Stopping rule
The adaptive rule certifies the sign of the *noisy* mean decision function (k_noisy expectation), not the exact classifier. Either (i) apply the correction (f-b)/a + b inside the estimator (a treated as fixed, adds its uncertainty ~ sd(k_id)/0.75 ~ 0.002 relative), or (ii) accept that the reference is the noisy classifier. Disagreement with exact is bounded to points with |f_exact| < ~0.014 -- the same near-tie set the rule cannot resolve anyway.
## Limitations
FakeFez snapshot calibration (not today's device); Aer noise lacks drift, crosstalk, leakage, TLS; layout [0,1] chosen as first coupling edge, not best qubits; opt-level 3 on a parametrized template (not per-point specialised). Calibration pairs come from training data (exact k classical). Literature (from memory, not fetched; no web used): Hubregtsen et al., arXiv:2105.02276 (kernel noise / mitigation); Thanasilp et al., arXiv:2208.11060 (noise-induced concentration).
## Reproduce
cd workstreams/w4_noise; P=~/.venvs/qhack-qml/bin/python; $P run_noise.py (~1 min, writes raw.npz); $P analyze.py > analyze_out.txt; $P extra.py. Functions: mitigation.py.
