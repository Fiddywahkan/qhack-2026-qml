# W3 kernel circuit report (FakeFez, dev data only, reps=2 unless noted)

Compiled cost (t1.py; 20 random dev x SV pairs, physical pair (0,1), seed_transpiler=7; duration = qc.estimate_duration(target), includes measurement; level 2 and 3 produce identical gate counts):

| variant | CZ | 1q (rz/sx/x) | depth | est. duration | dev acc | #SV | oracle mean shots/circuit | mean abs(k_noisy-k_exact), best pair (4,5) / median pair (40,41) |
|---|---|---|---|---|---|---|---|---|
| original L0 | 8 | 84 | 77 | 2.71 us | 0.9951 | 20 | 38.6 | 0.0109 / 0.0233 |
| original L1 | 5.6 | 40.1 | 35.8 | 2.30 us | same | 20 | 38.6 | 0.0091 / 0.0215 |
| original L2/L3 | 2.75 | 26.6 | 18.8 | 1.97 us | same | 20 | 38.6 | L3: 0.0070 / 0.0191 |
| merged L0 | 6 | 69 | 61 | 2.45 us | same | 20 | 38.6 | 0.0094 / 0.0217 |
| merged L1 | 5.6 | 39.9 | 35.6 | 2.30 us | same | 20 | 38.6 | not run |
| merged L2/L3 | 2.75 | 26.6 | 18.8 | 1.97 us | same | 20 | 38.6 | L3: 0.0070 / 0.0191 |

Findings
1. Transpiler: levels 0 and 1 do NOT fully merge the middle layers (L0 original 8 CZ; L1 5.6). Levels 2/3 resynthesize the whole 2-qubit unitary (KAK), giving <=3 CZ (mean 2.75, compact) for any variant, so hand-merging gives no benefit at L2/L3. Hand-merged circuit (fuse D(x_i)^dag D(x) into one RZZ + phases) saves 2 CZ at L0 only. Equivalence check (t1.py): max abs k difference 1.1e-15 over 90 pairs, reps 1-3.
2. Reps trade-off (t2.py, dev n=813, accuracy of sign(f) vs true labels):
   - reps=1: acc 0.8708, 26 SV, oracle 174.5 shots/circ (4536/pred), 30 near-ties
   - reps=2: acc 0.9951, 20 SV, 38.6 shots/circ (772/pred), 5 near-ties
   - reps=3: acc 0.7454, 37 SV, 374.9 shots/circ (13871/pred), 81 near-ties
   reps=2 is best on every axis (accuracy, SV count, shots); keep it. Since any 2-qubit circuit compiles to <=3 CZ, deeper reps cost little extra gate depth at L3 but hurt the classifier.
3. Layout (t3.py, 176 FakeFez edges, merged L3, 3 sample circuits; estimated error = 1-prod(1-gate err incl. readout) + idle/duration decoherence term 0.5*dur*(1/T1+1/T2) per qubit): best pair (4,5) 0.042 (gate 0.0205 + deco 0.0215); next (1,2) 0.051, (142,143) 0.051; median pair (40,41) 0.0915 (gate 0.0652 + deco 0.0263). Best is ~2.2x lower. The pair choice must be redone on the live ibm_fez calibration (FakeFez snapshot is stale).
4. Noise preview (t4.py; AerSimulator.from_backend(FakeFez), 20 random pairs, 4096 shots, seeds 100+n): see table. Shot-noise floor of this measurement is ~0.0047 (t5.py), so excess noise at best pair, L3 is ~0.002-0.003 vs ~0.014 at median pair; transpile level matters (L0 0.0109 -> L3 0.0070 at best pair); bias is ~0 (-0.0005) at best pair, +0.0015 at median pair. Preview is only a depolarizing/readout/thermal Aer model, not real drift/crosstalk.

Recommendation: merged circuit (or original, identical at L3), optimization_level=3, seed_transpiler=7, pair (4,5) per FakeFez (re-select on live calibration with t3.py logic), reps=2. Function: workstreams/w3_circuit/recommended.py (kernel_circuit, built on circuits.py).

Limitations: estimated duration is from the fake backend target (includes measurement); layout ranking used only 3 sample circuits and an approximate error model; the worst-pair score >1 reflects a bad-qubit decoherence term and is not a probability; noise-run not done for original L2 or merged L1; single transpiler seed. A circuit emitting P(00)=k via classically precomputed k would be 0-CZ but defeats the quantum kernel, so not considered.

Reproduce (from this dir): ~/.venvs/qhack-qml/bin/python t1.py; t2.py; t3.py (writes layout.pkl); t4.py (needs layout.pkl, ~8 min); t5.py. (t34.py is the superseded combined script.)
