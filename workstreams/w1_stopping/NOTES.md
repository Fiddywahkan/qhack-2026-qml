# W1: stopping rule (shots per circuit), dev data only

All numbers: `final.py` (sim.py engine, calib.py boundaries), seed 0, schedule 8..2048 doubling per circuit, L=9 checks, plug-in
sigma_hat = sqrt(sum w_i^2 kt_i(1-kt_i)/s_i), kt=(c+1)/(s+2), stop when |f_hat| >= z_j*sigma_hat, forced decision at cap.
Primary table is T=1000 trials (`table_T1000.md`, json alongside); T=200 (harness default) in `table_T200.md`. 813 dev points, 5 near-ties.
Worst-disagreement at T=200 is noisy (e.g. Bonf d=.3 gives .055 at T=200 and .058 at T=1000; Pocock d=.2 .060/.066; Bonf d=.2 .030/.040).

Cost columns are per prediction (20 circuits). "kind": principled = threshold set in advance from delta=0.05 only; tuned = delta (or z) chosen by looking at dev.

| rule | params | kind | mean shots | p95 | mean rounds | worst non-near | strict | agree | dev acc |
|---|---|---|---|---|---|---|---|---|---|
| v2 Bonferroni (reproduced) | d=.05, z=2.539 | principled | 1093 | 2560 | 2.40 | 0.006 | yes | 0.9984 | 0.9957 |
| Bonferroni | d=.1 | tuned | 942 | 2560 | 2.19 | 0.014 | yes | 0.9982 | 0.9956 |
| Bonferroni | d=.2 | tuned | 787 | 2560 | 1.98 | 0.040 | yes | 0.9978 | 0.9953 |
| Bonferroni | d=.3 | tuned | 689 | 1280 | 1.85 | 0.058 | NO | 0.9973 | 0.9950 |
| Pocock-calibrated const z | d=.05, z=2.419 | principled | 1020 | 2560 | 2.30 | 0.018 | yes | 0.9983 | 0.9957 |
| Pocock-calibrated const z | d=.1, z=2.119 | tuned | 844 | 2560 | 2.06 | 0.035 | yes | 0.9980 | 0.9955 |
| Pocock-calibrated const z | d=.2, z=1.762 | tuned | 654 | 1280 | 1.80 | 0.066 | NO | 0.9972 | 0.9948 |
| Lan-DeMets Pocock-type spend | d=.05 | principled | 1316 | 5120 | 2.96 | 0.007 | yes | 0.9985 | 0.9957 |
| Lan-DeMets OBF-type spend | d=.05 | principled | 2852 | 10240 | 4.22 | 0.004 | yes | 0.9984 | 0.9958 |
| poly spend delta_j ~ j^-2 | d=.05 | principled | 1218 | 5120 | 2.18 | 0.026 | yes | 0.9978 | 0.9952 |
| Robbins normal-mixture CS | d=.05, rho=1 | principled | 1700 | 5120 | 3.04 | 0.001 | yes | 0.9985 | 0.9957 |
| v2 Bonferroni + Neyman | d=.05 | principled | 976 | 2570 | 2.28 | 0.007 | yes | 0.9984 | 0.9957 |
| Bonferroni + Neyman | d=.1 | tuned | 847 | 2570 | 2.10 | 0.012 | yes | 0.9983 | 0.9956 |
| Bonferroni + Neyman | d=.2 | tuned | 712 | 1292 | 1.90 | 0.038 | yes | 0.9979 | 0.9953 |
| **Pocock-calibrated + Neyman** | d=.05, z=2.419 | principled | **912** | 2570 | 2.19 | 0.011 | yes | 0.9984 | 0.9957 |
| Pocock-calibrated + Neyman | d=.1, z=2.119 | tuned | 767 | 2568 | 1.98 | 0.027 | yes | 0.9981 | 0.9954 |
| Pocock-calibrated + Neyman | d=.2 | tuned | 597 | 1291 | 1.74 | 0.066 | NO | 0.9972 | 0.9948 |
| Lan-DeMets Pocock-type + Neyman | d=.05 | principled | 1176 | 2572 | 2.80 | 0.002 | yes | 0.9986 | 0.9958 |

Rounds distribution (fraction of predictions stopping at round 1..9) is in the last column of `table_T1000.md`. For the recommended rule (T=200): 39/31/15/8/4/2/1/0/1 %.
Cap (2048/circuit): only the near-tie points hit it (4 of 5 near-ties in >50% of trials for most rules; every non-near point stops before cap in nearly all trials). Non-near points are never "uncertain" in practice.

## Findings
- Principled beats principled: Pocock-type group-sequential calibration (constant z, chosen so that the one-sided crossing prob. under the Gaussian null with information 1,2,4,..,256 is delta) gives z=2.419 vs Bonferroni 2.539: -7% shots (1020 vs 1093), fewer rounds, same strict result. With Neyman: 912 (-17% vs v2, 1093).
- Tuned Bonferroni wins on raw cost but not on principle: d=.2 -> 787 (712 with Neyman) and still meets strict, but d=.3 fails (worst .058). The principled Pocock+Neyman (912) does not beat tuned Bonf d=.1/.2 + Neyman (847/712). So answer to the key question: NO, no principled rule here beats tuned-Bonferroni on cost; the margin is the slack between nominal delta and realised error (per-point error ~1% at nominal 5%, because the strict criterion is on realised wrong-sign rate and plug-in boundary crossing at true f != 0 is far below the null-crossing bound).
- Time-uniform rules (normal-mixture, OBF/Pocock spending, poly spending) are more expensive here: with only 9 pre-specified checks and a doubling schedule, discrete group-sequential is tighter than a continuous-time boundary. Back-loaded spending (OBF) is very costly (2852) because early looks are nearly blind.
- Finer schedule / SPRT: not run (SPRT for composite H0:f=0 is covered by the mixture CS above; no gain). Empirical-Bernstein not run; plug-in variance used throughout.

## Recommendation
Principled: stop when |f_hat| >= 2.419*sigma_hat (Pocock constant boundary calibrated for K=9 doubling looks, one-sided delta=0.05; Jennison & Turnbull, Group Sequential Methods, 2000, ch.2; Pocock 1977), with Neyman allocation after the pilot round. ~912 shots/pred, 2.2 rounds. If tuned-on-dev is acceptable: Bonferroni+Neyman delta=0.1-0.2 -> 712-847, but delta=.2 has worst .038 (margin thin; .3 fails), and test-split generalisation is unverified.

## Limitations
- sigma_hat is a plug-in estimate, not a true variance; at 8 shots kt can be misestimated, so delta is nominal, not guaranteed. Validity here is empirical (T=1000, worst over 808 points).
- z=2.419 was computed by Monte Carlo (4e6 paths, seed 123) of a Gaussian random walk, not exact integration (+-~0.01).
- Strict criterion with ~808 points x T trials is a max over noisy rates; the .05 line has MC noise (~+-.01 at T=1000).
- Noiseless simulation only; hardware noise biases k (not modelled). Each round = separate job on hardware (mean ~2.2).
- Neyman: cumulative per-circuit targets via ceil of budget share, using kt from previous rounds; first round uniform 8.

## Reproduce
```
cd workstreams/w1_stopping
~/.venvs/qhack-qml/bin/python final.py 1000   # ~3 min, prints table
~/.venvs/qhack-qml/bin/python final.py 200
~/.venvs/qhack-qml/bin/python calib.py        # boundaries
```
