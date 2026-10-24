"""Cost of the distribution-free (Hoeffding + union bound) sign certificate for direct measurement vs the empirical
Pocock rule, noiseless, dev only. Noise-aware sign certificate (eps=0): stop when
|b + mean| > R*sqrt(log(2L/delta)/(2n)), R = lam_max - lam_min; abstain at the cap."""
import numpy as np
d = np.load("direct_probs.npz"); lam, P = d["lam"], d["P_dev"]
import sys; sys.path.insert(0, "../harness"); import qhack_harness as h
dd = h.load(); b, f, near = dd["b"], dd["f_dev"], dd["near"]
T, delta, sched = 200, 0.05, 16 * 2 ** np.arange(9); L = len(sched); R = lam.max() - lam.min()
rng = np.random.default_rng(0)
cnt = np.zeros((T, len(P), 4), np.int64); active = np.ones((T, len(P)), bool)
f_stop = np.zeros((T, len(P))); s_stop = np.full((T, len(P)), sched[-1]); abstain = np.zeros((T, len(P)), bool)
for l, n in enumerate(sched):
    cnt += np.stack([rng.multinomial(n - (sched[l - 1] if l else 0), p, size=T) for p in P], axis=1)
    fh = b + (cnt / n) @ lam
    stop = active & (np.abs(fh) > R * np.sqrt(np.log(2 * L / delta) / (2 * n)))
    f_stop[stop] = fh[stop]; s_stop[stop] = n; active &= ~stop
abstain = active
wrong = (np.sign(f_stop) != np.sign(f)[None]) & ~abstain
print(f"R = {R:.3f}; certified rule on dev (noiseless, T={T}): shots/pred mean {s_stop.mean():.0f} "
      f"(non-near {s_stop[:, ~near].mean():.0f}); abstain rate {abstain.mean():.4f} "
      f"(non-near {abstain[:, ~near].mean():.4f}); wrong certified signs {wrong.sum()} of {(~abstain).sum()}")
