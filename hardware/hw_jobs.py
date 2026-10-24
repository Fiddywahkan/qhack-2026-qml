"""Hardware jobs for the adaptive-shots QML experiment (IBM Quantum, Heron).

Nothing here submits to a real device unless called with --submit. Without it, the same job code runs locally
on the FakeFez simulator (a snapshot noise model of ibm_fez), so every job can be dry-run first.

  python hw_jobs.py plan                                  # print circuit/shot counts for every job
  python hw_jobs.py timing [--submit --backend ibm_fez]   # job set 1: QPU-time probe
  python hw_jobs.py main --pool 2048 [--submit ...]       # (superseded) collect one shot pool for replay
  python hw_jobs.py live [--rep-delay 50e-6] [--submit ...] # job set 2: live fixed-1024, fixed-48, live adaptive

Results (job ids, per-shot bits, usage seconds) are saved under hardware/runs/.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from qiskit.circuit.library import zz_feature_map
from qiskit_machine_learning.datasets import ad_hoc_data
from qiskit_machine_learning.kernels import FidelityStatevectorKernel
from qiskit_machine_learning.utils import algorithm_globals

DATA_SEED, SPLIT_SEED, TRANSPILE_SEED, POINT_SEED = 12345, 2026, 7, 31
RUNS = Path(__file__).parent / "runs"


# ---------- data and model (identical to the notebook, sections 1-2) ----------
def _matches(A, B, tol=1e-9):
    return (np.abs(A[:, None, :] - B[None, :, :]).max(axis=2) < tol).any(axis=1)


def load_model():
    algorithm_globals.random_seed = DATA_SEED
    X_train, y_train, _, _ = ad_hoc_data(training_size=20, test_size=10, n=2, gap=0.3, one_hot=False)
    algorithm_globals.random_seed = DATA_SEED
    X_big, y_big, _, _ = ad_hoc_data(training_size=1000, test_size=0, n=2, gap=0.3, one_hot=False)
    keep = ~_matches(X_big, X_train)
    X_rest, y_rest = X_big[keep], y_big[keep]
    _, first = np.unique(X_rest.round(9), axis=0, return_index=True)
    first = np.sort(first)
    X_dev, X_test, y_dev, y_test = train_test_split(
        X_rest[first], y_rest[first], test_size=0.5, stratify=y_rest[first], random_state=SPLIT_SEED)
    qk = FidelityStatevectorKernel(feature_map=zz_feature_map(feature_dimension=2, reps=2))
    svm = SVC(kernel="precomputed").fit(qk.evaluate(X_train), y_train)
    sv = svm.support_
    return dict(X_dev=X_dev, X_test=X_test, y_test=y_test, X_sv=X_train[sv], w=svm.dual_coef_[0],
                b=svm.intercept_[0], qk=qk, classes=svm.classes_)


def choose_points(m, n_random=30, n_hard=10):
    """Test points for the hardware run. Selection uses features and the trained model only (|f| from the
    exact kernel), never labels: n_random uniformly at random + the n_hard closest to the boundary that are
    still decidable within the cap (oracle s* <= 2048 at z=2.5)."""
    K = m["qk"].evaluate(m["X_test"], m["X_sv"])
    f = K @ m["w"] + m["b"]
    s_star = 2.5**2 * ((m["w"]**2) * K * (1 - K)).sum(1) / f**2
    rng = np.random.default_rng(POINT_SEED)
    rand = rng.choice(len(f), n_random, replace=False)
    order = [i for i in np.argsort(np.abs(f)) if s_star[i] <= 2048 and i not in set(rand)]
    idx = np.concatenate([rand, order[:n_hard]])
    return idx, np.r_[np.zeros(n_random, bool), np.ones(n_hard, bool)]


# ---------- circuits ----------
def merged_kernel(x, xi, reps=2):
    """U(xi)^dag U(x) for zz_feature_map(2, reps) with the middle diagonal layers fused (identical k; verified
    to 1.1e-15 by workstream 3). Measures both qubits; k = P('00')."""
    def diag(v):
        return 2 * v[0], 2 * v[1], 2 * (np.pi - v[0]) * (np.pi - v[1])
    qc = QuantumCircuit(2)
    for _ in range(reps - 1):
        qc.h([0, 1]); a, b, t = diag(x); qc.p(a, 0); qc.p(b, 1); qc.rzz(t, 0, 1)
    qc.h([0, 1])
    (a, b, t), (a2, b2, t2) = diag(x), diag(xi)
    qc.p(a - a2, 0); qc.p(b - b2, 1); qc.rzz(t - t2, 0, 1)
    qc.h([0, 1])
    for _ in range(reps - 1):
        a, b, t = diag(xi); qc.rzz(-t, 0, 1); qc.p(-a, 0); qc.p(-b, 1); qc.h([0, 1])
    qc.measure_all()
    return qc


def direct_observable(m):
    """Decision observable M = sum_i w_i |psi(x_i)><psi(x_i)| and its eigenbasis (Nakayama et al., arXiv 2408.03000;
    variance proof in the vault, qml-literature/direct-observable-proof.md). f(x) = b + <psi(x)|M|psi(x)>."""
    from qiskit.quantum_info import Statevector
    fmap = zz_feature_map(feature_dimension=2, reps=2)
    psi = lambda x: Statevector(fmap.assign_parameters(x)).data
    M = sum(wi * np.outer(psi(xi), psi(xi).conj()) for wi, xi in zip(m["w"], m["X_sv"]))
    lam, V = np.linalg.eigh(M)
    return lam, V, fmap


def build_direct(X, backend, pair, V, fmap):
    from qiskit.circuit.library import UnitaryGate
    out = []
    for x in X:
        qc = QuantumCircuit(2)
        qc.compose(fmap.assign_parameters(x), inplace=True)
        qc.append(UnitaryGate(V.conj().T, label="Vdag"), [0, 1])
        qc.measure_all()
        out.append(transpile(qc, backend, initial_layout=list(pair), optimization_level=3,
                             seed_transpiler=TRANSPILE_SEED))
    return out


def best_pair(backend):
    """Lowest (2-qubit gate error + both readout errors) edge in the backend's current calibration."""
    tgt = backend.target
    op = "cz" if "cz" in tgt.operation_names else "ecr"
    best = None
    for (q0, q1), props in tgt[op].items():
        if props is None or props.error is None:
            continue
        ro = sum(tgt["measure"][(q,)].error or 0 for q in (q0, q1))
        score = props.error + ro
        if best is None or score < best[0]:
            best = (score, (q0, q1), props.error, ro)
    return best


def build(pairs_xy, backend, pair):
    return [transpile(merged_kernel(x, xi), backend, initial_layout=list(pair),
                      optimization_level=3, seed_transpiler=TRANSPILE_SEED) for x, xi in pairs_xy]


# ---------- job sets ----------
def timing_pubs(m):
    """QPU-time probe on dev points (no test data): fits usage = overhead + per_circuit*C + per_shot*shots."""
    x0 = m["X_dev"][0]
    small = [(x0, xi) for xi in m["X_sv"]]                                  # 20 circuits
    large = [(x, xi) for x in m["X_dev"][1:11] for xi in m["X_sv"]]          # 200 circuits
    # last job: same as t20x2000 but with a 50 us repetition delay instead of the 250 us default, to test whether
    # the per-shot wait (not the ~2 us circuit) dominates QPU time, and what it costs in kernel accuracy
    return [("t20x100", small, 100, None), ("t20x2000", small, 2000, None), ("t200x100", large, 100, None),
            ("t20x2000_rd50", small, 2000, 50e-6)]


def main_pubs(m, n_random=30, n_hard=10):
    idx, hard = choose_points(m, n_random, n_hard)
    kernel = [(m["X_test"][i], xi) for i in idx for xi in m["X_sv"]]
    return idx, hard, kernel, calibration_pairs(m)


def calibration_pairs(m, n=5):
    """Noise-calibration circuits: support-vector pairs whose exact k is known classically, spread over the k range.
    (Self-pairs k(x, x) = 1 are useless after optimization_level=3: the compiler cancels U^dag U to an empty
    circuit, so they would see readout error only.) Known exact k is available here because 2 qubits are
    classically simulable; on a classically hard kernel, structure-matched mirror circuits would be needed."""
    K = m["qk"].evaluate(m["X_sv"])
    iu = np.triu_indices(len(K), 1)
    order = np.argsort(K[iu])
    picks = order[np.linspace(0, len(order) - 1, n + 2).astype(int)[1:-1]]
    return [(m["X_sv"][iu[0][j]], m["X_sv"][iu[1][j]]) for j in picks]


def get_backend(args):
    if args.submit:
        from qiskit_ibm_runtime import QiskitRuntimeService
        service = QiskitRuntimeService()                                   # uses the saved account
        if args.backend == "least_busy":
            return service.least_busy(operational=True, simulator=False)
        return service.backend(args.backend)
    from qiskit_ibm_runtime.fake_provider import FakeFez
    return FakeFez()


def run_job(backend, circuits, shots, tag, args, rep_delay=None):
    from qiskit_ibm_runtime import SamplerV2
    sampler = SamplerV2(mode=backend)
    if args.submit:
        sampler.options.max_execution_time = args.max_seconds           # IBM stops the job past this
        if rep_delay is not None:
            sampler.options.execution.rep_delay = rep_delay
    t0 = time.time()
    job = sampler.run([(c,) for c in circuits], shots=shots)
    print(f"[{tag}] submitted job {job.job_id()} ({len(circuits)} circuits x {shots} shots)", flush=True)
    res = job.result()
    bits = np.stack([r.data.meas.array[:, 0] for r in res])               # packed 2-bit outcomes per shot
    usage = None
    if args.submit:
        try:
            usage = job.usage()
        except Exception as e:                                             # usage can lag the result
            print("usage not yet available:", e)
    return dict(tag=tag, job_id=job.job_id(), circuits=len(circuits), shots=shots, rep_delay=rep_delay, usage_s=usage,
                wall_s=time.time() - t0), bits


def run_live(backend, m, pair, args):
    """Job set 2. Same test points, three policies, every one executed for real so QPU usage is measured:
      fixed1024  baseline for the strict target (+ 5 calibration circuits to measure the noise scale a afterwards)
      fixed48    baseline for the lenient target
      adaptive   Pocock z=2.419 + Neyman allocation, one job per round with all still-undecided points batched,
                 capped at 1024 shots/circuit (8 looks) so it never spends more per point than fixed1024;
                 runs with --rep-delay if given (the reset-wait lever), baselines use the device default."""
    from pipeline import adaptive_predict
    from qiskit_ibm_runtime import SamplerV2
    idx, hard, kern, cal = main_pubs(m, args.n_random, args.n_hard)
    n = len(idx)
    circs = build(kern + cal, backend, pair)
    n_sv = len(m["X_sv"])
    grid = np.empty((n, n_sv), dtype=object)        # explicit fill: np.array() would unpack each circuit's gates
    for i in range(n):
        for j in range(n_sv):
            grid[i, j] = circs[i * n_sv + j]
    meta = dict(backend=backend.name, real=args.submit, pair=pair, rep_delay_adaptive=args.rep_delay,
                arms=args.arms, test_idx=idx.tolist(), hard=hard.tolist(), jobs=[],
                cz=[c.count_ops().get("cz", 0) for c in circs])
    arrays = dict(test_idx=idx, hard=hard)

    arms = args.arms.split(",")
    if "fixed1024" in arms:
        info, bits = run_job(backend, circs, 1024, "fixed1024", args)
        meta["jobs"].append(info); arrays["fixed1024_bits"] = bits
    if "fixed48" in arms:
        info, bits = run_job(backend, circs[:len(kern)], 48, "fixed48", args)
        meta["jobs"].append(info); arrays["fixed48_bits"] = bits

    for j in meta["jobs"]:
        j.update(arm=j["tag"], total_shots=j["shots"] * j["circuits"])

    def sampler_for_round(rep_delay):
        sampler = SamplerV2(mode=backend)
        if args.submit:
            sampler.options.max_execution_time = args.max_seconds
            if rep_delay is not None:
                sampler.options.execution.rep_delay = rep_delay
        return sampler

    def log_round(arm, job, n_circ, total_shots, t0, rep_delay):
        usage = None
        if args.submit:
            try:
                usage = job.usage()
            except Exception as e:
                print("usage not yet available:", e)
        k = len([j for j in meta["jobs"] if j.get("arm") == arm]) + 1
        meta["jobs"].append(dict(tag=f"{arm}_round{k}", arm=arm, job_id=job.job_id(), circuits=n_circ,
                                 shots=None, total_shots=total_shots, rep_delay=rep_delay, usage_s=usage,
                                 wall_s=time.time() - t0))
        print(f"[{arm} round {k}] job {job.job_id()} ({n_circ} circuits, {total_shots} shots)", flush=True)

    def live_counts(add):
        out = np.zeros(add.shape, np.int64)
        todo = np.argwhere(add > 0)
        if len(todo) == 0:
            return out
        t0 = time.time()
        job = sampler_for_round(None).run([(grid[i, j], None, int(add[i, j])) for i, j in todo])
        res = job.result()
        for (i, j), r in zip(todo, res):
            out[i, j] = int((r.data.meas.array[:, 0] == 0).sum())
        log_round("adaptive", job, len(todo), int(add.sum()), t0, None)
        return out

    if "adaptive" in arms:
        res = adaptive_predict(live_counts, m["w"], m["b"], a=1.0, z=2.419, s0=8, looks=8, neyman=True,
                               batch_shape=(n,))
        arrays.update(adaptive_f=res["f"], adaptive_shots=res["shots"], adaptive_rounds=res["rounds"],
                      adaptive_capped=res["capped"], adaptive_s=res["s"])

    # Direct decision-observable measurement (one circuit per point). Arms isolate one lever each:
    #   direct_d2_rd250  doubling 16..4096 (9 looks), default 250 us reset wait
    #   direct_d2_rd50   same, 50 us reset wait            -> reset-wait lever
    #   direct_d4_rd50   x4 growth 16..4096 (5 looks)      -> fewer rounds (job overhead) lever
    #   direct_fixed2048_rd50  one job, 2048 shots each    -> adaptive vs a single fixed job
    # z = 2.419 was calibrated for 9 doubling looks; with 5 looks it is conservative.
    direct_arms = {"direct_d2_rd250": (16 * 2 ** np.arange(9), None),
                   "direct_d2_rd50": (16 * 2 ** np.arange(9), 50e-6),
                   "direct_d4_rd50": (16 * 4 ** np.arange(5), 50e-6),
                   "direct_fixed2048_rd50": (np.array([2048]), 50e-6)}
    if any(a in arms for a in direct_arms):
        lam, V, fmap = direct_observable(m)
        dcircs = build_direct(m["X_test"][idx], backend, pair, V, fmap)
        meta["direct_cz"] = [c.count_ops().get("cz", 0) for c in dcircs]
        arrays.update(direct_lam=lam, direct_V=V)
    for arm, (sched, rd) in direct_arms.items():
        if arm not in arms:
            continue
        fixed = len(sched) == 1
        cnt = np.zeros((n, 4), np.int64); active = np.ones(n, bool)
        f_stop = np.zeros(n); s_stop = np.zeros(n, int); capped = np.zeros(n, bool)
        round_counts = []
        for l, tot in enumerate(sched):
            add = tot - (sched[l - 1] if l else 0)
            todo = np.flatnonzero(active)
            t0 = time.time()
            job = sampler_for_round(rd).run([(dcircs[i], None, int(add)) for i in todo])
            res = job.result()
            rc = np.zeros((n, 4), np.int64)
            for i, r in zip(todo, res):
                rc[i] = np.bincount(r.data.meas.array[:, 0], minlength=4)[:4]
            cnt += rc; round_counts.append(rc)
            log_round(arm, job, len(todo), int(add * len(todo)), t0, rd)
            pt = (cnt + 0.5) / (tot + 2.0)
            f = m["b"] + (cnt / tot) @ lam
            sig = np.sqrt((pt @ lam**2 - (pt @ lam) ** 2) / tot)
            last = l == len(sched) - 1
            decided = np.abs(f) >= 2.419 * sig
            stop = active & (decided | last)
            capped |= active & last & ~decided & ~np.bool_(fixed)
            f_stop[stop] = f[stop]; s_stop[stop] = tot; active &= ~stop
            if not active.any():
                break
        arrays.update({f"{arm}_f": f_stop, f"{arm}_shots": s_stop, f"{arm}_capped": capped,
                       f"{arm}_counts": cnt, f"{arm}_round_counts": np.stack(round_counts)})

    save("live", meta, arrays)
    for arm in dict.fromkeys(j["arm"] for j in meta["jobs"]):
        jobs = [j for j in meta["jobs"] if j["arm"] == arm]
        print(f"{arm:24s} jobs {len(jobs)}  total shots {sum(j['total_shots'] for j in jobs)}")


def usage_report(path):
    """Re-query billed usage for every job id in a saved run (usage can lag job completion)."""
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = QiskitRuntimeService()
    meta = json.loads(Path(path).read_text())
    for j in meta["jobs"]:
        j["usage_s"] = service.job(j["job_id"]).usage()
        print(f"{j['tag']:28s} {j.get('total_shots', j['shots']):>8} shots  usage {j['usage_s']} s")
    Path(path).write_text(json.dumps(meta, indent=2, default=str))


def save(name, meta, arrays):
    out = RUNS if meta.get("real") else RUNS / "dryrun"          # dry runs (FakeFez) never mix with real data
    out.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    np.savez_compressed(out / f"{name}-{stamp}.npz", **arrays)
    (out / f"{name}-{stamp}.json").write_text(json.dumps(meta, indent=2, default=str))
    print("saved", out / f"{name}-{stamp}.*")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["plan", "timing", "main", "live", "usage"])
    ap.add_argument("--run", help="for 'usage': path to a saved run .json")
    ap.add_argument("--n-random", type=int, default=15)
    ap.add_argument("--n-hard", type=int, default=5)
    ap.add_argument("--rep-delay", type=float, default=None, help="seconds; adaptive/direct arms only (e.g. 50e-6)")
    ap.add_argument("--arms", default="fixed1024,fixed48,adaptive,direct_d2_rd250,direct_d2_rd50,direct_d4_rd50,"
                    "direct_fixed2048_rd50", help="comma list of live arms to run")
    ap.add_argument("--submit", action="store_true", help="run on real IBM hardware (needs saved account)")
    ap.add_argument("--backend", default="ibm_fez", help='device name, or "least_busy"')
    ap.add_argument("--max-seconds", type=int, default=60, help="hard QPU-time cap per job")
    ap.add_argument("--pool", type=int, default=2048, help="shots per circuit in the main collection job")
    args = ap.parse_args()
    if args.what == "usage":
        return usage_report(args.run)
    m = load_model()
    if args.what == "plan":
        for tag, xy, shots, rd in timing_pubs(m):
            print(f"timing {tag}: {len(xy)} circuits x {shots} shots = {len(xy) * shots} shots"
                  + (f", rep_delay {rd * 1e6:.0f} us" if rd else ""))
        idx, hard, kern, cal = main_pubs(m, args.n_random, args.n_hard)
        n_c = len(kern) + len(cal)
        print(f"live: {len(idx)} test points; fixed1024 {n_c} circuits x 1024 = {n_c * 1024} shots; "
              f"fixed48 {len(kern)} x 48 = {len(kern) * 48} shots; adaptive: rounds of batched jobs, "
              f"cap 1024/circuit (worst case {len(kern) * 1024} shots)")
        print(f"main: {len(idx)} test points ({hard.sum()} near-boundary) x 20 SVs = {len(kern)} kernel circuits"
              f" + {len(cal)} calibration circuits, x {args.pool} shots = {(len(kern) + len(cal)) * args.pool} shots")
        return
    backend = get_backend(args)
    score, pair, gate_err, ro_err = best_pair(backend)
    print(f"backend {backend.name}: best pair {pair} (2q err {gate_err:.4f}, readout sum {ro_err:.4f})")
    meta = dict(backend=backend.name, real=args.submit, pair=pair, gate_err=gate_err, ro_err=ro_err,
                transpile_seed=TRANSPILE_SEED, max_exec_s=args.max_seconds, jobs=[])
    if args.what == "live":
        return run_live(backend, m, pair, args)
    if args.what == "timing":
        arrays = {}
        for tag, xy, shots, rd in timing_pubs(m):
            info, bits = run_job(backend, build(xy, backend, pair), shots, tag, args, rep_delay=rd)
            meta["jobs"].append(info); arrays[tag] = bits
        save("timing", meta, arrays)
    else:
        idx, hard, kern, cal = main_pubs(m, args.n_random, args.n_hard)
        circs = build(kern + cal, backend, pair)
        meta["cz_counts"] = [c.count_ops().get("cz", 0) for c in circs]
        info, bits = run_job(backend, circs, args.pool, "main", args)
        meta["jobs"].append(info)
        meta.update(test_idx=idx.tolist(), hard=hard.tolist(), n_kernel=len(kern), n_calibration=len(cal))
        save("main", meta, dict(bits=bits, test_idx=idx, hard=hard))


if __name__ == "__main__":
    sys.exit(main())
