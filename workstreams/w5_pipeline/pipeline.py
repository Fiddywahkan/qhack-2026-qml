"""Reusable adaptive prediction: Pocock-boundary sequential test + Neyman allocation + affine-noise rescale.
adaptive_predict(counts_fn, w, b, a, ...) works on a BATCH of predictions (leading axes) so it can drive
simulation (counts_fn = Binomial draws) or replay of real hardware shot pools (counts_fn = draw from a pool).

counts_fn(alloc): alloc int array (..., m) = extra shots requested per SV circuit for each prediction (0 for
  predictions that already stopped); must return int array (..., m) of '00' counts out of exactly alloc shots.
Mitigation: noisy k = a*k + c; sum(w)=0 so c cancels and f_mit = (sum_i w_i khat_i)/a + b. The standard error
  is computed on the same (mitigated) scale: sigma_mit = sqrt(sum w_i^2 kt(1-kt)/s_i)/a, kt=(c+1)/(s+2).
  Stop when |f_mit| >= z*sigma_mit. a=1 reproduces the unmitigated rule (certifies the noisy classifier).
Returns dict: f (decision value at stop), shots (total per prediction), rounds, capped (bool), s (per-SV shots)."""
import numpy as np


def adaptive_predict(counts_fn, w, b, a=1.0, z=2.419, s0=8, looks=9, neyman=True, batch_shape=()):
    w = np.asarray(w, float); m = len(w)
    shp = tuple(batch_shape) + (m,)
    s = np.zeros(shp, np.int64); c = np.zeros(shp, np.int64)
    active = np.ones(shp[:-1], bool)
    f_stop = np.zeros(shp[:-1]); rounds = np.zeros(shp[:-1], int); capped = np.zeros(shp[:-1], bool)
    for j in range(looks):
        B = m * s0 * 2**j                       # cumulative total budget at look j (per-circuit average s0*2^j)
        if j == 0 or not neyman:
            target = np.full(shp, s0 * 2**j, float)
        else:
            kt = (c + 1) / (s + 2); p = np.abs(w) * np.sqrt(kt * (1 - kt)) + 1e-12
            target = np.maximum(B * p / p.sum(-1, keepdims=True), s)
        add = np.where(active[..., None], np.maximum(np.ceil(target - s), 0), 0).astype(np.int64)
        c = c + counts_fn(add); s = s + add
        rounds += active
        kt = (c + 1) / (s + 2)
        f = ((c / np.maximum(s, 1)) @ w) / a + b
        sig = np.sqrt((w**2 * kt * (1 - kt) / np.maximum(s, 1)).sum(-1)) / a
        stop = active & (np.abs(f) >= z * sig)
        if j == looks - 1:
            capped = active & ~stop; stop = active.copy()
        f_stop[stop] = f[stop]; active &= ~stop
    return dict(f=f_stop, shots=s.sum(-1), rounds=rounds, capped=capped, s=s)


def estimate_a_identity(k_id_hat):
    """a from circuits with exact k=1 under full-depolarising-to-1/4 model."""
    return (np.mean(k_id_hat) - 0.25) / 0.75


def estimate_a_calpairs(k_noisy_hat, k_exact):
    """a from calibration pairs with known exact k (needed when identity circuits compile to nothing):
    k_noisy-1/4 = a (k-1/4), LS through the origin."""
    u = np.asarray(k_exact) - 0.25; v = np.asarray(k_noisy_hat) - 0.25
    return float((u * v).sum() / (u * u).sum())
