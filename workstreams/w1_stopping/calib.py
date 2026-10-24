# Group-sequential boundaries for Gaussian canonical joint distribution, information V_j = 2^j (doubling), one-sided.
# Z_j = S_j/sqrt(V_j), S = Brownian motion at V_j. Calibrated by Monte Carlo (seed fixed). Ref: Jennison & Turnbull 2000;
# Lan & DeMets 1983 (alpha spending).
import numpy as np
from scipy.stats import norm
L=9
def paths(n=4_000_000, seed=123):
    rng=np.random.default_rng(seed)
    V=2.0**np.arange(L); dV=np.diff(np.r_[0,V])
    S=np.cumsum(rng.standard_normal((n,L))*np.sqrt(dV),axis=1)
    return S/np.sqrt(V)
Z=paths()
def pocock(delta):
    lo,hi=1.0,4.0
    for _ in range(40):
        mid=(lo+hi)/2
        p=(Z>=mid).any(1).mean()
        lo,hi=(mid,hi) if p>delta else (lo,mid)
    return hi
def spend_bounds(delta,spend):  # spend(t)->cumulative alpha, t=V_j/V_L
    t=2.0**np.arange(L)/2.0**(L-1); cum=spend(t,delta); zs=[]
    alive=np.ones(len(Z),bool); prev=0
    for j in range(L):
        target=(cum[j]-prev)*len(Z)
        zz=np.sort(Z[alive,j])[::-1]
        k=int(round(target)); k=max(k,1)
        z=zz[k-1] if k<=len(zz) else zz[-1]
        zs.append(z); alive&=~(Z[:,j]>=z); prev=cum[j]
    return np.array(zs)
obf=lambda t,a: 2-2*norm.cdf(norm.isf(a/2)/np.sqrt(t))   # OBF-type (two-sided form -> one-sided alpha a: use a directly)
def obf1(t,a): return 2*(1-norm.cdf(norm.isf(a/2)/np.sqrt(t)))*(a/ (a))  # cumulative, ends at a
def pock1(t,a): return a*np.log(1+(np.e-1)*t)
if __name__=="__main__":
    for dl in (0.05,0.1,0.2): print('pocock',dl,round(pocock(dl),3),'bonf',round(norm.isf(dl/L),3))
    print('OBF-type d=.05',np.round(spend_bounds(0.05,lambda t,a:2*(1-norm.cdf(norm.isf(a/2)/np.sqrt(t)))),2))
    print('Pock-type d=.05',np.round(spend_bounds(0.05,pock1),2))
