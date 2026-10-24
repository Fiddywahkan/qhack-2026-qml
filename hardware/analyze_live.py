"""Analyze real-hardware live runs (hardware/runs/live-*.json/.npz with real=true), merging multiple invocations.

Per arm: predictions vs the exact classifier and vs true labels on the same 20 test points, total shots,
billed QPU seconds (IBM job usage), jobs/rounds, wall time. Also the real-chip kernel noise fit from fixed1024.
  python analyze_live.py            # prints tables, writes runs/live_summary.json and ../figures/06_hardware_*.png
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import hw_jobs as h

RUNS = Path(__file__).parent / "runs"
FIG = Path(__file__).parent.parent / "figures"


def load_real_runs():
    metas, arrays = [], {}
    for js in sorted(glob.glob(str(RUNS / "live-*.json"))):
        meta = json.loads(Path(js).read_text())
        if not meta.get("real"):
            continue
        metas.append(meta)
        z = np.load(js.replace(".json", ".npz"))
        arrays.update({k: z[k] for k in z.files})
    assert metas, "no real live runs found"
    idx = [tuple(mm["test_idx"]) for mm in metas]
    assert len(set(idx)) == 1, "runs used different test points"
    return metas, arrays


def main():
    metas, A = load_real_runs()
    m = h.load_model()
    idx = np.array(metas[0]["test_idx"]); hard = np.array(metas[0]["hard"])
    n, n_sv = len(idx), len(m["X_sv"])
    K = m["qk"].evaluate(m["X_test"][idx], m["X_sv"])
    f_ex = K @ m["w"] + m["b"]
    y = m["y_test"][idx]
    label = lambda f: np.where(f > 0, m["classes"][1], m["classes"][0])
    jobs = [j for mm in metas for j in mm["jobs"]]
    print(f"backend {metas[0]['backend']}, pair {metas[0]['pair']}, {n} test points ({hard.sum()} near-boundary)")
    print(f"exact classifier on these points: accuracy {np.mean(label(f_ex) == y):.3f}")

    rows = {}
    def add(arm, f_hat):
        js = [j for j in jobs if j.get("arm", j["tag"]) == arm]
        usage = [j["usage_s"] for j in js]
        rows[arm] = dict(agree=float(np.mean(np.sign(f_hat) == np.sign(f_ex))),
                         agree_random=float(np.mean((np.sign(f_hat) == np.sign(f_ex))[~hard])),
                         agree_hard=float(np.mean((np.sign(f_hat) == np.sign(f_ex))[hard])),
                         accuracy=float(np.mean(label(f_hat) == y)),
                         jobs=len(js), shots=int(sum(j["total_shots"] for j in js)),
                         usage_s=None if any(u is None for u in usage) else float(sum(usage)),
                         wall_s=float(sum(j["wall_s"] for j in js)))

    if "fixed1024_bits" in A:
        bits = A["fixed1024_bits"]
        k_hw = (bits[: n * n_sv] == 0).mean(1).reshape(n, n_sv)
        add("fixed1024", k_hw @ m["w"] + m["b"])
        slope, off = np.polyfit(K.ravel(), k_hw.ravel(), 1)
        cal = h.calibration_pairs(m)
        k_cal_ex = np.array([m["qk"].evaluate(a[None], b[None])[0, 0] for a, b in cal])
        k_cal_hw = (bits[n * n_sv:] == 0).mean(1)
        print(f"real-chip kernel noise (fixed1024, {K.size} circuits): k_hw = {slope:.3f} k + {off:.4f}; "
              f"mean |k_hw - k| {np.abs(k_hw - K).mean():.4f}; calibration pairs: {np.round(k_cal_hw, 3)} vs exact {np.round(k_cal_ex, 3)}")
        rows["fixed1024"]["noise_fit"] = [float(slope), float(off)]
    if "fixed48_bits" in A:
        k48 = (A["fixed48_bits"] == 0).mean(1).reshape(n, n_sv)
        add("fixed48", k48 @ m["w"] + m["b"])
    if "adaptive_f" in A:
        add("adaptive", A["adaptive_f"])
    for arm in ("direct_d2_rd250", "direct_d2_rd50", "direct_d4_rd50", "direct_fixed2048_rd50"):
        if f"{arm}_f" in A:
            add(arm, A[f"{arm}_f"])
            rows[arm]["capped"] = int(A[f"{arm}_capped"].sum())
    if "direct_fixed2048_rd50_counts" in A:                       # direct-measurement noise: 2048-shot estimate vs exact
        lam = A["direct_lam"]; c = A["direct_fixed2048_rd50_counts"]
        f_d = m["b"] + (c / c.sum(1, keepdims=True)) @ lam
        sl, of = np.polyfit(f_ex - m["b"], f_d - m["b"], 1)
        print(f"real-chip direct measurement (2048 shots, 50 us): f_hw - b = {sl:.3f}(f - b) + {of:.4f}; "
              f"mean |f_hw - f| {np.abs(f_d - f_ex).mean():.4f}")

    print(f"\n{'arm':24s} {'agree':>6s} {'(rand)':>6s} {'(hard)':>6s} {'acc':>6s} {'jobs':>5s} {'shots':>8s} "
          f"{'QPU s':>6s} {'QPU s/pred':>10s} {'shots/pred':>10s}")
    for arm, r in rows.items():
        u = r["usage_s"]
        print(f"{arm:24s} {r['agree']:6.2f} {r['agree_random']:6.2f} {r['agree_hard']:6.2f} {r['accuracy']:6.2f} "
              f"{r['jobs']:5d} {r['shots']:8d} {u if u is not None else float('nan'):6.0f} "
              f"{(u / n) if u is not None else float('nan'):10.2f} {r['shots'] / n:10.0f}")
    (RUNS / "live_summary.json").write_text(json.dumps(dict(rows=rows, n=n, backend=metas[0]["backend"],
                                                            pair=metas[0]["pair"]), indent=2))
    plot(rows, n)


def plot(rows, n):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = {"fixed1024": "fixed 1024, 20 circuits (strict baseline)", "fixed48": "fixed 48, 20 circuits (lenient baseline)",
             "adaptive": "adaptive, 20 circuits", "direct_d2_rd250": "direct adaptive, 9 rounds, 250 µs wait",
             "direct_d2_rd50": "direct adaptive, 9 rounds, 50 µs wait", "direct_d4_rd50": "direct adaptive, 5 rounds, 50 µs wait",
             "direct_fixed2048_rd50": "direct fixed 2048, 1 job, 50 µs wait"}
    arms = [a for a in names if a in rows and rows[a]["usage_s"] is not None]
    if not arms:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
    y = np.arange(len(arms))[::-1]
    for ax, key, xlabel in [(axes[0], "shots", "shots per prediction (log)"), (axes[1], "usage_s", "billed QPU seconds per prediction (log)")]:
        vals = [rows[a][key] / n for a in arms]
        ax.barh(y, vals, color=["0.55" if a.startswith("fixed") else "#2a78d6" for a in arms], height=0.6, edgecolor="white")
        for yi, v in zip(y, vals):
            ax.text(v * 1.12, yi, f"{v:,.2f}" if v < 10 else f"{v:,.0f}", va="center", fontsize=8.5, color="0.2")
        ax.set_xscale("log"); ax.set_xlabel(xlabel); ax.grid(axis="y", visible=False)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(y, [names[a] + f"  (agree {rows[a]['agree']:.0%})" for a in arms], fontsize=8.5)
    axes[1].set_xlim(right=max(rows[a]["usage_s"] / n for a in arms) * 4)
    axes[0].set_xlim(right=max(rows[a]["shots"] / n for a in arms) * 6)
    fig.suptitle(f"Real hardware ({n} test points): shots vs billed QPU time per prediction", fontsize=10)
    fig.tight_layout(); FIG.mkdir(exist_ok=True); fig.savefig(FIG / "06_hardware_cost.png", dpi=110)
    print("saved figures/06_hardware_cost.png")


if __name__ == "__main__":
    main()
