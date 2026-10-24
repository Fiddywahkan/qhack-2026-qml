import sys, numpy as np
from scipy.stats import norm
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness'))
from qhack_harness import load, evaluate
S0, CAPS = 8, 2048
L = 9  # checks at 8,16,...,2048 per circuit

def run(d, zfun, neyman=False, T=200, seed=0):
    rng = np.random.default_rng(seed)
    K = np.broadcast_to(d["K_dev"], (T,)+d["K_dev"].shape)  # T,N,m
    w, b = d["w"], d["b"]; m = len(w)
    shape = K.shape[:2]
    s = np.zeros(K.shape); c = np.zeros(K.shape)
    active = np.ones(shape, bool); fstop = np.zeros(shape); rounds = np.zeros(shape); capped = np.zeros(shape, bool)
    for j in range(L):
        B = m*S0*2**j   # cumulative total shot budget at check j
        if j == 0 or not neyman:
            target = np.full(K.shape, S0*2**j, float)
        else:
            kt = (c+1)/(s+2); p = np.abs(w)*np.sqrt(kt*(1-kt)) + 1e-12
            target = B*p/p.sum(-1, keepdims=True)
            target = np.maximum(target, s)
        add = np.where(active[...,None], np.ceil(target - s), 0).astype(int)
        add = np.maximum(add, 0)
        c += rng.binomial(add, K); s += add
        rounds += active
        kt = (c+1)/(s+2)
        f = (c/np.maximum(s,1))@w + b if False else ((c+0.5*0)/s)@w + b
        sig = np.sqrt((w**2*kt*(1-kt)/s).sum(-1))
        z = zfun(j, sig)
        stop = active & (np.abs(f) >= z*sig)
        if j == L-1:
            capped = active & ~stop; stop = active
        fstop[stop] = f[stop]; active &= ~stop
    tot = s.sum(-1)
    return fstop, tot, rounds, capped

def zbonf(delta, Lc=L):
    z = norm.isf(delta/Lc); return lambda j, sig: z
def zspend(deltas):
    zs = norm.isf(np.asarray(deltas)); return lambda j, sig: zs[j]
def zpoly(delta, p=2):  # delta_j = delta * j^-p / zeta(p) , j=1..inf
    from scipy.special import zeta
    return zspend([delta*(j+1)**-p/zeta(p) for j in range(L)])
def zmix(delta, rho):
    # Robbins normal-mixture boundary (two-sided): z^2 = (1+rho/V)(log(1+V/rho)+2log(2/delta)), V=1/sig^2 (information)
    def f(j, sig):
        V = 1/sig**2
        return np.sqrt((1+rho/V)*(np.log1p(V/rho)+2*np.log(2/delta)))
    return f

def row(name, params, kind, d, zf, ney=False, T=200, seed=0):
    f, tot, rd, cap = run(d, zf, ney, T, seed)
    e = evaluate(d, f, tot)
    caprate = cap.mean(0)
    return dict(rule=name, params=params, kind=kind, mean=e['mean_shots_per_prediction'], p95=e['p95_shots_per_prediction'],
        rounds=rd.mean(), worst=e['worst_far'], strict=e['meets_strict'], agree=e['agree'], acc=e['accuracy'],
        cap_pts=int((caprate>0.5).sum()), cap_frac=cap.mean(), near=int(d['near'].sum()), n=len(caprate),
        rounds_hist=np.bincount(rd.ravel().astype(int), minlength=L+1)[1:]/rd.size)
