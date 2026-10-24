# The cheapest trustworthy predictions from a quantum-kernel classifier

QHack — Interstate Quantum Hackathon for Beginners, Fall 2026 (MTSU) · QML track · Kieran McCrorey (UAH)

**Question.** A quantum-kernel classifier gets each prediction from many repeated circuit runs (shots), and shots cost QPU time.
What is the cheapest way to get predictions that still agree with the exact classifier, and does it hold on a real IBM quantum computer?

| Result | Baseline | Best |
|---|---|---|
| Shots per prediction, held-out test set (simulation, 814 points) | 20,480 | **65** (direct measurement + Pocock stopping), same accuracy |
| Billed QPU seconds per prediction, `ibm_fez` (20 test points, batch average) | 5.50 s (fixed 1024, strict baseline) | **0.25 s** (direct measurement, one job), all predictions agree with the exact classifier |

Read that second row with its control: in the same run `fixed48` used **0.40 s** and also reproduced all 20 predictions, so the direct arm is 1.6× faster than fixed 48 on this batch, not 22×. The case for 1024 shots as the *strict* baseline comes from the 814-point simulation, where fixed 48 leaves a boundary point wrong 33% of the time. Reset wait and job shape changed along with the measurement architecture, so the usage gain is not attributable to the direct construction alone.

On real hardware, the policy with the fewest shots was *not* the cheapest in QPU time: every job is billed at least ~2 s, so one batched job beat multi-round adaptive stopping.

## Contents

| Path | What |
|---|---|
| [`report.pdf`](report.pdf) | **The project report.** |
| `qhack_adaptive_shots.ipynb` | **The technical record.** Every method, what worked and what didn't, all results and figures. Starts with the install cell; runs top to bottom in ~1 min without IBM access (hardware results load from `hardware/runs/`). |
| `figures/` | Figures produced by the notebook (axes labelled). |
| `hardware/hw_jobs.py` | IBM hardware jobs (timing probe, live A/B). Dry-runs on the FakeFez simulator unless `--submit`. |
| `hardware/analyze_live.py`, `hardware/pipeline.py` | Hardware analysis; adaptive stopping + Neyman allocation. |
| `hardware/runs/` | Raw per-shot results, job ids and billed usage of the real runs (`dryrun/`: simulator rehearsals; the ones before 14:19 predate the circuit-grid fix in notebook §16). |
| `workstreams/` | Workstreams (W1 stopping rule, W2 circuits per prediction, W3 circuit/compilation, W4 noise, W5 integration, direct observable), each with `NOTES.md` (method and results) and scripts. |

## Reproduce

```
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute qhack_adaptive_shots.ipynb
```
