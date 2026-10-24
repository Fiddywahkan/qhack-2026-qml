import sys, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness')); sys.path.insert(0,'.')
from qhack_harness import load
from mitigation import *
d=load(); R=np.load('raw.npz'); sub=R['sub']; Pd=R['Pd']; A=R['A']
w,b=d['w'],d['b']; Kex=d['K_dev'][sub]; fex=d['f_dev'][sub]; y=d['y_dev'][sub]; near=d['near'][sub]; cl=d['classes']
raw=Pd[...,0]
ro=ro_mitigate_p00_vec(Pd,A)
def stats(name,extra,K,a=1,c=0):
    f=decision(K,w,b,a,c); e=np.abs(correct_kernel(K,a,c)-Kex)
    ag=(np.sign(f)==np.sign(fex)); pred=np.where(f>0,cl[1],cl[0])
    print(f"| {name} | {extra} | {e.mean():.4f} | {e.max():.4f} | {ag.mean():.3f} (far {ag[~near].mean():.3f}) | {(pred==y).mean():.3f} |")
print("exact classifier dev acc on subset:", (np.where(fex>0,cl[1],cl[0])==y).mean(), "n=",len(sub))
# bias structure
for nm,K in [('raw',raw),('ro',ro)]:
    a_,c_=np.polyfit(Kex.ravel(),K.ravel(),1); res=K.ravel()-(a_*Kex.ravel()+c_)
    print(nm,'fit k_noisy=a*k+c: a=%.4f c=%.4f resid std=%.4f max=%.4f; corr=%.5f'%(a_,c_,res.std(),np.abs(res).max(),np.corrcoef(Kex.ravel(),K.ravel())[0,1]))
    q=np.polyfit(Kex.ravel(),K.ravel(),2); print('  quad coeffs',q)
print('shot sd approx', np.sqrt(0.25/20000))
Pc=R['Pc']; kc=R['kc_exact']; Pid=R['Pid']
rc=ro_mitigate_p00_vec(Pc,A); rid=ro_mitigate_p00_vec(Pid,A)
print('identity raw k:',Pid[:,0],' ro-mit:',rid)
print('| setting | extra circuits (shots each) | mean|dk| | max|dk| | agree w/ exact | dev acc |')
stats('exact (ref)','-',Kex) if False else None
stats('unmitigated','0',raw)
stats('RO mitigation','4 (cal)',ro)
a1,c1=identity_affine(Pid[:,0]); stats('depol rescale, identity only (no RO)','4 identity',raw,a1,c1)
a2,c2=identity_affine(rid); stats('RO + depol rescale (identity)','4 cal + 4 id',ro,a2,c2)
a3,c3=fit_affine(Pc[:,0],kc); stats('depol affine fit on 40 cal pairs (no RO)','40',raw,a3,c3)
a4,c4=fit_affine(rc,kc); stats('RO + affine fit (40 cal pairs)','4 + 40',ro,a4,c4)
print('params: id-only a,c',a1,c1,' RO+id',a2,c2,' fit',a3,c3,' RO+fit',a4,c4)
# few-cal robustness: fit with n pairs
rng=np.random.default_rng(0)
for n in (5,10,20):
    ac=[]
    for _ in range(200):
        i=rng.choice(40,n,replace=False); a_,c_=fit_affine(rc[i],kc[i])
        f=decision(ro,w,b,a_,c_); ac.append((np.sign(f)==np.sign(fex)).mean())
    print(f'RO+affine fit with {n} cal pairs: agreement mean {np.mean(ac):.3f} min {np.min(ac):.3f}')
# only sign-threshold shift: f_noisy sign vs threshold. unmitigated bias in f
fn=raw@w+b; print('f bias: mean(fn-fex)=%.4f ; sum w=%.4f b=%.4f'%((fn-fex).mean(),w.sum(),b))
print('disagree unmitig: n=',(np.sign(fn)!=np.sign(fex)).sum(),' |fex| of disagreers',np.sort(np.abs(fex[np.sign(fn)!=np.sign(fex)]))[:10])
