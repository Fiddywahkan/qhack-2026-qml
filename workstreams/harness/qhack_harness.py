"""Shared evaluation harness for the QHack adaptive-shots experiments (DEV SET ONLY).

Reproduces the notebook's sections 1-3 exactly (same seeds) and exposes:
  d = load()  -> dict with X_train, y_train, X_dev, y_dev, X_sv, sv, w, b, K_dev (dev x SV exact kernel),
                 f_dev (exact decision values), near (bool, near-tie dev points: oracle s*(z=2.5) > 2048),
                 svm, qkernel, fmap, classes
  evaluate(d, f_stop, shots_per_prediction) -> metrics dict
     f_stop: (trials, n_dev) estimated decision values at the moment the policy stops
     shots_per_prediction: (trials, n_dev) TOTAL shots spent on that prediction (all circuits, all rounds)

Strict criterion (the project's main target): no non-near-tie dev point disagrees with the exact
classifier in more than 5% of trials. Baselines on this harness (T=200 trials, 20 SVs):
  fixed shots:             1024 shots/circuit  = 20480 shots/prediction
  adaptive v2 (Bonferroni, start 8, delta=0.05, uniform):  ~54.5 shots/circuit = ~1090 shots/prediction
  adaptive v3 (+ Neyman allocation):                       ~49   shots/circuit = ~985  shots/prediction
  oracle (exact k, optimal allocation, z=2.5):             ~34   shots/circuit

THE TEST SPLIT IS DELIBERATELY NOT EXPOSED. Do not regenerate or use it.
"""
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from qiskit.circuit.library import zz_feature_map
from qiskit_machine_learning.datasets import ad_hoc_data
from qiskit_machine_learning.utils import algorithm_globals
from qiskit_machine_learning.kernels import FidelityStatevectorKernel

DATA_SEED, SPLIT_SEED, CAP, STRICT = 12345, 2026, 2048, 0.05


def _matches(A, B, tol=1e-9):
    return (np.abs(A[:, None, :] - B[None, :, :]).max(axis=2) < tol).any(axis=1)


def load(reps=2):
    algorithm_globals.random_seed = DATA_SEED
    X_train, y_train, _, _ = ad_hoc_data(training_size=20, test_size=10, n=2, gap=0.3, one_hot=False)
    algorithm_globals.random_seed = DATA_SEED
    X_big, y_big, _, _ = ad_hoc_data(training_size=1000, test_size=0, n=2, gap=0.3, one_hot=False)
    keep = ~_matches(X_big, X_train)
    X_rest, y_rest = X_big[keep], y_big[keep]
    _, first = np.unique(X_rest.round(9), axis=0, return_index=True)
    first = np.sort(first)
    X_pool, y_pool = X_rest[first], y_rest[first]
    X_dev, _X_test, y_dev, _y_test = train_test_split(
        X_pool, y_pool, test_size=0.5, stratify=y_pool, random_state=SPLIT_SEED)
    del _X_test, _y_test
    fmap = zz_feature_map(feature_dimension=2, reps=reps)
    qkernel = FidelityStatevectorKernel(feature_map=fmap)
    svm = SVC(kernel="precomputed").fit(qkernel.evaluate(X_train), y_train)
    sv = svm.support_
    w, b = svm.dual_coef_[0], svm.intercept_[0]
    K_dev = qkernel.evaluate(X_dev, X_train[sv])
    f_dev = K_dev @ w + b
    s_star = 2.5**2 * ((w**2) * K_dev * (1 - K_dev)).sum(1) / f_dev**2
    return dict(X_train=X_train, y_train=y_train, X_dev=X_dev, y_dev=y_dev, X_sv=X_train[sv], sv=sv,
                w=w, b=b, K_dev=K_dev, f_dev=f_dev, near=s_star > CAP, svm=svm, qkernel=qkernel,
                fmap=fmap, classes=svm.classes_)


def evaluate(d, f_stop, shots_per_prediction):
    near = d["near"]
    dis = (np.sign(f_stop) != np.sign(d["f_dev"])[None]).mean(0)
    pred = np.where(f_stop > 0, d["classes"][1], d["classes"][0])
    return dict(agree=1 - dis.mean(), agree_far=1 - dis[~near].mean(), worst_far=dis[~near].max(),
                n_far_over_strict=int((dis[~near] > STRICT).sum()),
                meets_strict=bool(dis[~near].max() <= STRICT),
                accuracy=(pred == d["y_dev"][None]).mean(),
                mean_shots_per_prediction=float(np.mean(shots_per_prediction)),
                p95_shots_per_prediction=float(np.percentile(shots_per_prediction, 95)))


if __name__ == "__main__":
    d = load()
    T = 200
    rng = np.random.default_rng(0)
    for s in (512, 1024):
        F = (rng.binomial(s, np.broadcast_to(d["K_dev"], (T, *d["K_dev"].shape))) / s) @ d["w"] + d["b"]
        print(s, evaluate(d, F, np.full(F.shape, s * len(d["sv"]))))
