# Factor 2: circuits per prediction (dev only, T=200, seed-fixed)
Reproduce: `cd <this dir> && ~/.venvs/qhack-qml/bin/python main.py > out.txt` (~2 min; writes rows.json with every number; common.py holds the simulator).
All rows use v2 rule: shots 8,8,16,32,...(cumulative doubling), z=2.54 (Bonferroni over L=9). 813 dev pts, 5 near-ties.
"worst" = worst non-near-tie disagreement vs the *method's own exact classifier* (strict as the harness defines it). "agree orig / acc" = vs ORIGINAL exact f and vs true labels.

| method | params | circ | shots | rounds | worst (own ref) | strict | worst vs ORIG exact | agree ORIG | dev acc |
|---|---|---|---|---|---|---|---|---|---|
| v2 reproduced | 20 SV | 20 | 1190 | 2.68 | 0.005 | yes | 0.005 | 0.9985 | 0.9957 |
| lazy, |w|-order | a0=4,x2 | 19.6 | 1592 | 4.1 | 0 | yes | 0 | 0.9985 | 0.9959 |
| lazy | a0=8,x2 | 19.9 | 1309 | 3.5 | 0.01 | yes | 0.01 | 0.9983 | 0.9958 |
| lazy | a0=12,x1.5 | 19.5 | 1266 | 3.2 | 0 | yes | 0 | 0.9986 | 0.9957 |
| lazy | a0=12,x2 | 20.0 | 1209 | 3.0 | 0 | yes | 0 | 0.9985 | 0.9958 |
| top-16 |w|, refit b | m=16 | 16 | 813 | 2.8 | 0.005 | yes | 0.995 | 0.9951 | 1.000 |
| top-16, LSQ refit to f_orig on train | m=16 | 16 | 1113 | 2.9 | 0.01 | yes | 0.02 | 0.9982 | 0.9958 |
| top-12 LSQ refit | m=12 | 12 | 2731 | 4.5 | 0.025 | NO | 0.385 | 0.9954 | 0.9909 |
| top-12 refit b | m=12 | 12 | 2136 | 3.6 | 0.01 | yes | 1.0 | 0.951 | 0.956 |
| retrain C=0.5 | 26 SV | 26 | 1385 | 2.2 | 0.01 | yes | 0.985 | 0.9967 | 0.9924 |
| retrain C=3 | 15 SV | 15 | 1668 | 3.7 | 0.01 | yes | 0.57 | 0.996 | 0.9967 |
| retrain C=10 | 13 SV | 13 | 1879 | 4.1 | 0.015 | yes | 1.0 | 0.988 | 0.9907 |
Full table (m=6..16, C=0.1..10, more lazy settings) in out.txt/rows.json.
(Note: with T=200, 0.005 = 1 disagreement; the strict flag is within MC noise at 0.005-0.01, so "NO" at 0.025 is the only clear fail.)

## Findings
- **Lazy SV measurement: no gain (negative result).** Sorted |w| = 1.0 x7, then 0.94 ... 0.25: weights are flat (box constraint C=1 active), sum|w|=15.4. The deterministic bound on unmeasured SVs is ~sum|w_unmeas|, far wider than |f|, so no point can be decided before nearly all SVs are measured; lazy costs extra rounds and ~equal or more shots (1209-2476 vs 1190), and still ~20 circuits. Would only help with a heavy-tailed weight distribution (e.g. an L1 / large-C model).
- **Sparser classifier is the only lever that cuts circuits.** Best: drop 4 SVs (smallest |w|) and refit only b: 16 circuits, 813 shots (-32% shots, -20% circuits), dev acc 1.000 vs 0.9957 (3 pts of 813; m chosen on dev, so selection bias; do not read as improvement), agree with original 0.9951 (it is a *different* classifier: ~4 dev pts flip; worst_vs_orig 0.995 = one point that always differs). If the strict criterion must be vs the ORIGINAL classifier, this fails; vs its own exact classifier it passes. LSQ-refit variant (reduced-set, Burges 1996 style) keeps agreement 0.998 with original and passes vs both (worst 0.02 vs orig is 4/200, i.e. marginal fail vs orig), but costs 1113 shots (~v2) with 16 circuits. Below m=16 things degrade fast (12 circuits needs 2100-2700 shots, accuracy drops or strict fails). Retrained C=3 gives 15 SV, acc 0.9967, 1668 shots.
- Idea 3 (classical geometry prior) not run: the ad-hoc kernel is not distance-based, so a classical proxy is unlikely to identify relevant SVs and flat weights make ordering moot; no time.

## Error budget
Every row: per-check CI at level 0.05/9 (union over the 9 checks => total error <= 5% for the estimator-based sign); lazy rows add deterministic (zero-error) bounds for unmeasured SVs, so the budget is unchanged. Plug-in k~=(c+1)/(s+2) variance is heuristic (as in v2), not a guarantee; validity is checked empirically via the strict test.

## Winner
For circuits: top-16 + LSQ refit (16 circuits, ~1113 shots, passes against both references within noise) or top-16 refit-b (16 circuits, 813 shots, but a different classifier). Neither is a decisive win; the 20->16 circuit saving is modest. Decision: keep all 20 (the exact classifier is the reference). Caveats: single seed per row, T=200, m chosen on dev. Test split never touched.
