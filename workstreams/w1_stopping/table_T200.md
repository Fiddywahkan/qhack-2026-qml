| rule | params | kind | mean shots/pred | p95 | mean rounds | worst non-near disagr. | strict | agree | dev acc | cap pts(>50%) | rounds dist 1..9 |
| v2 Bonferroni | delta=.05,L=9 | principled | 1101 | 2560 | 2.40 | 0.015 | True | 0.9983 | 0.9957 | 5/5 near | [0.35, 0.26, 0.19, 0.1, 0.05, 0.02, 0.01, 0.0, 0.01] |
| Bonferroni | delta=.1 | tuned | 944 | 2560 | 2.20 | 0.025 | True | 0.9982 | 0.9955 | 4/5 near | [0.43, 0.25, 0.16, 0.08, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Bonferroni | delta=.2 | tuned | 791 | 2560 | 1.98 | 0.030 | True | 0.9978 | 0.9952 | 4/5 near | [0.51, 0.23, 0.13, 0.07, 0.03, 0.01, 0.01, 0.0, 0.01] |
| Bonferroni | delta=.3 | tuned | 684 | 1280 | 1.85 | 0.055 | False | 0.9975 | 0.9949 | 4/5 near | [0.57, 0.21, 0.12, 0.06, 0.03, 0.01, 0.0, 0.0, 0.0] |
| Pocock GS (MC-calibrated const z) | delta=0.05,z=2.419 | principled | 1017 | 2560 | 2.30 | 0.015 | True | 0.9984 | 0.9956 | 5/5 near | [0.39, 0.26, 0.17, 0.09, 0.05, 0.02, 0.01, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) | delta=0.1,z=2.119 | tuned | 843 | 2560 | 2.06 | 0.020 | True | 0.9981 | 0.9955 | 4/5 near | [0.48, 0.24, 0.15, 0.07, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) | delta=0.2,z=1.762 | tuned | 650 | 1280 | 1.80 | 0.060 | False | 0.9970 | 0.9949 | 4/5 near | [0.59, 0.2, 0.11, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0] |
| Lan-DeMets Pocock-type spend | delta=.05 | principled | 1315 | 5120 | 2.96 | 0.005 | True | 0.9985 | 0.9957 | 4/5 near | [0.15, 0.26, 0.29, 0.17, 0.08, 0.03, 0.01, 0.0, 0.01] |
| Lan-DeMets Pocock-type spend | delta=.1 | tuned | 1171 | 2560 | 2.80 | 0.010 | True | 0.9983 | 0.9959 | 4/5 near | [0.19, 0.27, 0.27, 0.15, 0.07, 0.03, 0.01, 0.0, 0.01] |
| Lan-DeMets OBF-type spend | delta=.05 | principled | 2858 | 10240 | 4.22 | 0.005 | True | 0.9984 | 0.9959 | 4/5 near | [0.02, 0.09, 0.24, 0.27, 0.16, 0.12, 0.07, 0.01, 0.01] |
| Lan-DeMets OBF-type spend | delta=.1 | tuned | 2661 | 10240 | 4.20 | 0.000 | True | 0.9985 | 0.9958 | 4/5 near | [0.02, 0.09, 0.24, 0.27, 0.16, 0.14, 0.05, 0.01, 0.01] |
| poly spend p=2 | delta=.05 | principled | 1216 | 5120 | 2.18 | 0.015 | True | 0.9981 | 0.9952 | 6/5 near | [0.55, 0.13, 0.12, 0.09, 0.05, 0.03, 0.01, 0.01, 0.01] |
| Robbins normal-mixture CS | delta=.05,rho=1 | principled | 1693 | 5120 | 3.04 | 0.005 | True | 0.9985 | 0.9958 | 6/5 near | [0.15, 0.27, 0.25, 0.16, 0.08, 0.04, 0.02, 0.01, 0.01] |
| v2 Bonferroni + Neyman | delta=.05,L=9 | principled | 978 | 2570 | 2.28 | 0.010 | True | 0.9985 | 0.9957 | 4/5 near | [0.35, 0.32, 0.16, 0.09, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Bonferroni + Neyman | delta=.1 | tuned | 848 | 2570 | 2.09 | 0.020 | True | 0.9985 | 0.9956 | 4/5 near | [0.43, 0.3, 0.14, 0.07, 0.03, 0.01, 0.01, 0.0, 0.01] |
| Bonferroni + Neyman | delta=.2 | tuned | 718 | 1293 | 1.90 | 0.030 | True | 0.9982 | 0.9954 | 4/5 near | [0.51, 0.27, 0.11, 0.06, 0.03, 0.01, 0.0, 0.0, 0.0] |
| Bonferroni + Neyman | delta=.3 | tuned | 628 | 1291 | 1.78 | 0.055 | False | 0.9974 | 0.9951 | 3/5 near | [0.57, 0.25, 0.1, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0] |
| Pocock GS (MC-calibrated const z) + Neyman | delta=0.05,z=2.419 | principled | 918 | 2570 | 2.19 | 0.010 | True | 0.9985 | 0.9955 | 4/5 near | [0.39, 0.31, 0.15, 0.08, 0.04, 0.02, 0.01, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) + Neyman | delta=0.1,z=2.119 | tuned | 766 | 2568 | 1.98 | 0.025 | True | 0.9981 | 0.9958 | 4/5 near | [0.48, 0.28, 0.12, 0.06, 0.03, 0.01, 0.0, 0.0, 0.01] |
| Pocock GS (MC-calibrated const z) + Neyman | delta=0.2,z=1.762 | tuned | 597 | 1291 | 1.74 | 0.060 | False | 0.9975 | 0.9951 | 2/5 near | [0.59, 0.24, 0.09, 0.05, 0.02, 0.01, 0.0, 0.0, 0.0] |
| Lan-DeMets Pocock-type spend + Neyman | delta=.05 | principled | 1172 | 2572 | 2.80 | 0.005 | True | 0.9986 | 0.9957 | 4/5 near | [0.15, 0.34, 0.26, 0.15, 0.07, 0.03, 0.01, 0.0, 0.01] |
| Lan-DeMets Pocock-type spend + Neyman | delta=.1 | tuned | 1057 | 2571 | 2.66 | 0.005 | True | 0.9985 | 0.9957 | 4/5 near | [0.19, 0.34, 0.25, 0.13, 0.06, 0.02, 0.01, 0.0, 0.01] |
