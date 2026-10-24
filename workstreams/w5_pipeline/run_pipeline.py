"""Pipeline ablation in noisy simulation. Needs probs.npz (from probs.py). T=200 trials, seed 0 stream per row."""
import sys, numpy as np
S=str(__import__('pathlib').Path(__file__).resolve().parents[1])
sys.path[:0]=[S+'/harness', '.']
from qhack_harness import load, evaluate
from pipeline import adaptive_predict, estimate_a_identity, estimate_a_calpairs
d=load(); R=np.load('probs.npz'); sub=R['sub']; w,b=d['w'],d['b']
ds=dict(d); ds['near']=d['near'][sub]; ds['f_dev']=d['f_dev'][sub]; ds['y_dev']=d['y_dev'][sub]
K=d['K_dev'][sub]; N=len(sub); T=200; FIX=1024; CALSH=1024
print('subset',N,'near',int(ds['near'].sum()),'sum(w)=%.2e'%w.sum())
def run_row(name, P, circuit, adaptive, mitig, seed, noisy_id=None):
    rng=np.random.default_rng(seed)
    # per-trial a estimate (shot noise included)
    if mitig:
        if circuit=='t':   # 4 identity circuits x 1024 shots
            kid=rng.binomial(CALSH, R['t_id'][None].repeat(T,0))/CALSH; a=((kid.mean(1)-.25)/.75)[:,None]; calc=4*CALSH
        else:              # 8 calibration pairs x 1024 shots, known exact k
            kc=rng.binomial(CALSH, R['w_cal'][None].repeat(T,0))/CALSH
            u=R['kcal']-.25; a=(((kc-.25)*u).sum(1)/(u*u).sum())[:,None]; calc=8*CALSH
    else: a=np.ones((T,1)); calc=0
    cf=lambda add: rng.binomial(add, np.broadcast_to(P,(T,N,20)))
    if adaptive:
        r=adaptive_predict(cf,w,b,a=a,batch_shape=(T,N)); f=r['f']; shots=r['shots']; rounds=r['rounds'].mean(); cap=r['capped'].mean()
    else:
        c=rng.binomial(FIX,np.broadcast_to(P,(T,N,20))); f=(c/FIX@w)/a+b; shots=np.full((T,N),FIX*20); rounds=1.0; cap=0.0
    e=evaluate(ds,f,shots)
    return dict(name=name,agree=e['agree'],agree_far=e['agree_far'],worst=e['worst_far'],strict=e['meets_strict'],nover=e['n_far_over_strict'],
        acc=e['accuracy'],mean=e['mean_shots_per_prediction'],p95=e['p95_shots_per_prediction'],rounds=rounds,cap=cap,
        dk=np.abs(P-K).mean(),a=(a.mean(),a.std()) if mitig else None,calc=calc,mfar=shots[:,~ds['near']].mean(),p95far=np.percentile(shots[:,~ds['near']],95),mnear=shots[:,ds['near']].mean())
rows=[
 run_row('BASELINE: template[0,1], fixed 1024, no mitig', R['t_P'],'t',False,False,1),
 run_row('OPTIMIZED: W3(4,5)+Pocock/Neyman+mitig', R['w_P'],'w',True,True,2),
 run_row('abl -W3 circuit (template[0,1]) +adaptive+mitig', R['t_P'],'t',True,True,3),
 run_row('abl -adaptive (fixed 1024) W3+mitig', R['w_P'],'w',False,True,4),
 run_row('abl -mitigation: W3+adaptive, a=1', R['w_P'],'w',True,False,5),
 run_row('extra: W3 only (fixed 1024, no mitig)', R['w_P'],'w',False,False,6),
 run_row('extra: template+adaptive, no mitig', R['t_P'],'t',True,False,7),
 run_row('NOISELESS adaptive (exact K)', K,'w',True,False,8),
]
print('| row | agree all | agree non-near | worst non-near dis | strict | #nonnear>5%% | dev acc | mean shots | p95 shots | rounds | capped frac | mean|dk| | a_hat (mean,sd) | calib shots/batch | non-near mean/p95 shots | near-tie mean shots |')
for r in rows:
    print('| %s | %.4f | %.4f | %.3f | %s | %d | %.4f | %.0f | %.0f | %.2f | %.3f | %.4f | %s | %d | %.0f / %.0f | %.0f |'%(r['name'],r['agree'],r['agree_far'],r['worst'],r['strict'],r['nover'],r['acc'],r['mean'],r['p95'],r['rounds'],r['cap'],r['dk'],
      'n/a' if r['a'] is None else '%.4f, %.4f'%r['a'],r['calc'],r['mfar'],r['p95far'],r['mnear']))
# true a (LS slope through 1/4) of each circuit family
for k,nm in (('t_P','template'),('w_P','W3')):
    u=(K-.25).ravel(); v=(R[k]-.25).ravel(); print(nm,'true a (LS thru 1/4) %.4f'%((u*v).sum()/(u*u).sum()),' identity-circuit k:',R[k[0]+'_id'].round(4),'cz id:',R[k[0]+'_cz'][N*20:N*20+4],'mean cz dev %.2f'%R[k[0]+'_cz'][:N*20].mean())
# near-tie / worst points detail for OPTIMIZED vs BASELINE: which non-near point is worst
