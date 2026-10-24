| rule | params | kind | mean shots/pred | p95 | mean rounds | worst non-near disagr. | strict | agree | dev acc | cap pts(>50%) | rounds dist 1..9 |
| v2 Bonferroni | delta=.05,L=9 | principled | 1093 | 2560 | 2.40 | 0.006 | True | 0.9984 | 0.9957 | 5/5 near | [0.35, 0.26, 0.19, 0.1, 0.05, 0.02, 0.01, 0.0, 0.01] |
| Bonferroni | delta=.1 | tuned | 942 | 2560 | 2.19 | 0.014 | True | 0.9982 | 0.9956 | 4/5 near | [0.43, 0.25, 0.16, 0.08, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Bonferroni | delta=.2 | tuned | 787 | 2560 | 1.98 | 0.040 | True | 0.9978 | 0.9953 | 4/5 near | [0.51, 0.23, 0.13, 0.07, 0.03, 0.01, 0.01, 0.0, 0.01] |
| Bonferroni | delta=.3 | tuned | 689 | 1280 | 1.85 | 0.058 | False | 0.9973 | 0.9950 | 4/5 near | [0.57, 0.21, 0.12, 0.06, 0.03, 0.01, 0.0, 0.0, 0.0] |
| Pocock GS (MC-calibrated const z) | delta=0.05,z=2.419 | principled | 1020 | 2560 | 2.30 | 0.018 | True | 0.9983 | 0.9957 | 5/5 near | [0.39, 0.26, 0.18, 0.09, 0.05, 0.02, 0.01, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) | delta=0.1,z=2.119 | tuned | 844 | 2560 | 2.06 | 0.035 | True | 0.9980 | 0.9955 | 4/5 near | [0.48, 0.24, 0.15, 0.07, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) | delta=0.2,z=1.762 | tuned | 654 | 1280 | 1.80 | 0.066 | False | 0.9972 | 0.9948 | 4/5 near | [0.59, 0.2, 0.11, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0] |
| Lan-DeMets Pocock-type spend | delta=.05 | principled | 1316 | 5120 | 2.96 | 0.007 | True | 0.9985 | 0.9957 | 4/5 near | [0.15, 0.26, 0.29, 0.17, 0.08, 0.03, 0.01, 0.0, 0.01] |
| Lan-DeMets Pocock-type spend | delta=.1 | tuned | 1175 | 2560 | 2.80 | 0.002 | True | 0.9984 | 0.9957 | 4/5 near | [0.19, 0.27, 0.27, 0.15, 0.07, 0.03, 0.01, 0.0, 0.01] |
| Lan-DeMets OBF-type spend | delta=.05 | principled | 2852 | 10240 | 4.22 | 0.004 | True | 0.9984 | 0.9958 | 4/5 near | [0.02, 0.1, 0.24, 0.27, 0.16, 0.12, 0.07, 0.01, 0.01] |
| Lan-DeMets OBF-type spend | delta=.1 | tuned | 2657 | 10240 | 4.19 | 0.003 | True | 0.9985 | 0.9957 | 4/5 near | [0.02, 0.1, 0.24, 0.27, 0.16, 0.14, 0.05, 0.01, 0.01] |
| poly spend p=2 | delta=.05 | principled | 1218 | 5120 | 2.18 | 0.026 | True | 0.9978 | 0.9952 | 6/5 near | [0.55, 0.13, 0.12, 0.09, 0.05, 0.03, 0.01, 0.01, 0.01] |
| Robbins normal-mixture CS | delta=.05,rho=1 | principled | 1700 | 5120 | 3.04 | 0.001 | True | 0.9985 | 0.9957 | 6/5 near | [0.16, 0.28, 0.25, 0.16, 0.08, 0.04, 0.02, 0.01, 0.01] |
| v2 Bonferroni + Neyman | delta=.05,L=9 | principled | 976 | 2570 | 2.28 | 0.007 | True | 0.9984 | 0.9957 | 4/5 near | [0.35, 0.32, 0.16, 0.09, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Bonferroni + Neyman | delta=.1 | tuned | 847 | 2570 | 2.10 | 0.012 | True | 0.9983 | 0.9956 | 4/5 near | [0.43, 0.3, 0.14, 0.07, 0.03, 0.01, 0.01, 0.0, 0.01] |
| Bonferroni + Neyman | delta=.2 | tuned | 712 | 1292 | 1.90 | 0.038 | True | 0.9979 | 0.9953 | 4/5 near | [0.51, 0.27, 0.11, 0.06, 0.03, 0.01, 0.0, 0.0, 0.0] |
| Bonferroni + Neyman | delta=.3 | tuned | 632 | 1291 | 1.78 | 0.057 | False | 0.9975 | 0.9950 | 4/5 near | [0.57, 0.25, 0.1, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0] |
| Pocock GS (MC-calibrated const z) + Neyman | delta=0.05,z=2.419 | principled | 912 | 2570 | 2.19 | 0.011 | True | 0.9984 | 0.9957 | 4/5 near | [0.39, 0.31, 0.15, 0.08, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) + Neyman | delta=0.1,z=2.119 | tuned | 767 | 2568 | 1.98 | 0.027 | True | 0.9981 | 0.9954 | 4/5 near | [0.48, 0.28, 0.12, 0.06, 0.03, 0.01, 0.0, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) + Neyman | delta=0.2,z=1.762 | tuned | 597 | 1291 | 1.74 | 0.066 | False | 0.9972 | 0.9948 | 3/5 near | [0.59, 0.24, 0.09, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0] |
| Lan-DeMets Pocock-type spend + Neyman | delta=.05 | principled | 1176 | 2572 | 2.80 | 0.002 | True | 0.9986 | 0.9958 | 4/5 near | [0.15, 0.34, 0.26, 0.14, 0.07, 0.03, 0.01, 0.0, 0.01] |
| Lan-DeMets Pocock-type spend + Neyman | delta=.1 | tuned | 1054 | 2571 | 2.66 | 0.004 | True | 0.9986 | 0.9958 | 4/5 near | [0.19, 0.35, 0.24, 0.13, 0.06, 0.02, 0.01, 0.0, 0.01] |
