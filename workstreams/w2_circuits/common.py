import sys
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness'))
import numpy as np
from scipy.stats import norm
from qhack_harness import load, evaluate, CAP
T=200; DELTA=0.05; S0=8; L=9

def run(w,b,K,order,a0,grow,seed,z=None,lazy=True,T=T,start=S0,rounds=L):
    """K: (n,m) exact kernel used only to draw binomials. order: (n,m) per-point SV order (indices) or (m,).
    Active set per round: first a_r SVs in order (a0, then *grow, capped m). Shots/round/active circuit: start*2^r (new draws,
    cumulative per SV). Stop when deterministic+CI interval excludes 0 (all L checks share delta -> Bonferroni)."""
    rng=np.random.default_rng(seed); n,m=K.shape
    z=z or norm.ppf(1-DELTA/rounds)
    order=np.broadcast_to(order,(n,m))
    wo=np.take_along_axis(np.broadcast_to(w,(n,m)),order,1); Ko=np.take_along_axis(K,order,1)
    c=np.zeros((T,n,m)); s=np.zeros((T,n,m))
    f=np.zeros((T,n)); shots=np.zeros((T,n)); circ=np.zeros((T,n)); rnd=np.zeros((T,n)); done=np.zeros((T,n),bool)
    a=a0
    for r in range(rounds):
        sr=start if r==0 else start*2**(r-1)
        act=np.arange(m)<a if lazy else np.ones(m,bool)
        idx=np.where(act)[0]
        # draw for active SVs of undone items
        cc=rng.binomial(sr,np.broadcast_to(Ko[:,idx],(T,n,len(idx))))
        live=~done
        c[:,:,idx]+=cc*live[...,None]; s[:,:,idx]+=sr*live[...,None]
        shots+=live*sr*len(idx); circ=np.maximum(circ, (s>0).sum(2)); rnd+=live
        meas=s>0
        kt=np.where(meas,(c+1)/(s+2),0.5)
        var=np.where(meas,wo**2*kt*(1-kt)/np.maximum(s,1),0).sum(2)
        zs=z*np.sqrt(var)
        fm=(wo*np.where(meas,kt,0)).sum(2)+b
        wu=np.where(meas,0,wo)
        lo=fm-zs+np.minimum(wu,0).sum(2); hi=fm+zs+np.maximum(wu,0).sum(2)
        fm_pt=fm+ (np.where(meas,0,wo)*0.5).sum(2)  # midpoint guess for unresolved
        dec=(lo>0)|(hi<0)
        if lazy: # when all SVs measured, deterministic part vanishes -> v2 rule
            pass
        new=dec&live
        f=np.where(new,np.where(lo>0,1.0,-1.0)*np.maximum(np.abs(fm),1e-12),f)
        last=(r==rounds-1)
        f=np.where(live&~new&last,fm,f)  # undecided at cap: sign of estimate (all measured? maybe not)
        done|=dec
        a=min(m,int(np.ceil(a*grow)))
    return f,shots,circ,rnd
