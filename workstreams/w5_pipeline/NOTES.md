# W5 pipeline ablation (FakeFez, Aer, dev only)
Method: noisy P(00) per circuit from `AerSimulator.from_backend(FakeFez(), method='density_matrix')`, **200000 shots/circuit** (not exact; sd ~1e-3 in k, same for all rows), 155 dev points (150 far, seed 0, + 5 near-ties) x 20 SVs, for (t) W4 parametrized template, layout [0,1], opt 3, seed_transpiler 11 (6 CZ) and (w) W3 circuit pair (4,5) (2.86 CZ mean). Then Binomial(s,p) shots, T=200 trials, seeds 1-8 per row. Per trial `a` is re-estimated from simulated calibration shots (1024 each).
Adaptive rule: `pipeline.py:adaptive_predict` (Pocock z=2.419, 9 looks, s0=8, Neyman, cap 2048/circuit avg; sigma and f on mitigated scale: f=(sum w khat)/a+b, sigma/a).
**Important finding:** W3's transpiled identity circuit k(x_sv,x_sv) compiles to **0 CZ, depth 1** (noise-free except readout; gives k=0.991 -> a=0.988, true a=0.974). So identity-circuit calibration is invalid for the W3 circuit. For W3 rows `a` is fit from 8 calibration circuits (train non-SV point vs SV with largest |k-1/4|, W3-transpiled, exact k known classically), k-1/4 = a(k_exact-1/4) LS. Template rows use the 4 identity circuits (6 CZ, valid).

Calibration shot cost is a fixed per-batch cost (4x1024 template, 8x1024 W3), NOT in shots/prediction.

Results: see run_out.txt (verbatim). Shots/prediction include the 5 near-ties (3% of this subset vs 0.6% of full dev), so use the non-near columns for the typical cost.

Reproduce: `cd w5_pipeline; P=~/.venvs/qhack-qml/bin/python; $P probs.py > probs.log` (~2 min, writes probs.npz); `$P run_pipeline.py > run_out.txt` (~4 s).
Limitations: FakeFez snapshot, Aer noise only (no drift/crosstalk), P(00) from 2e5 shots, one transpile seed, single 155-point subset, a shared across points within a trial, near-ties excluded from the strict criterion by definition.
