import sys,json
from common import *
from sklearn.svm import SVC
d=load(); w,b,K,fd=d["w"],d["b"],d["K_dev"],d["f_dev"]; m=len(w)
print("n_dev",K.shape, "n near",d["near"].sum(),"sum|w|",np.abs(w).sum(),"b",b)
rows=[]
def rep(name,params,f,sh,circ,rnd,dd=d,fo=None):
    e=evaluate(dd,f,sh)
    # agreement with ORIGINAL exact classifier
    ag=(np.sign(f)==np.sign(d["f_dev"])[None]).mean()
    dis_o=(np.sign(f)!=np.sign(d["f_dev"])[None]).mean(0)
    wo=dis_o[~d["near"]].max()
    r=dict(method=name,params=params,circ=float(circ.mean()),shots=e["mean_shots_per_prediction"],rounds=float(rnd.mean()),
     worst=float(e["worst_far"]),strict=e["meets_strict"],worst_vs_orig=float(wo),agree_orig=float(ag),acc=float(e["accuracy"]))
    rows.append(r); print(r,flush=True)
ordw=np.argsort(-np.abs(w))
# v2
f,sh,ci,rn=run(w,b,K,ordw,m,1,1,lazy=False); rep("v2 (all 20, doubling)","z=2.54",f,sh,ci,rn)
# lazy, |w| order
print('sorted |w|',np.round(np.abs(w)[ordw],2).tolist())
for a0 in (4,6,8,12,14,16,18):
  for g in ((1.5,2) if a0<14 else (2,)):
    f,sh,ci,rn=run(w,b,K,ordw,a0,g,2,lazy=True); rep("lazy |w|-order",f"a0={a0},grow={g}",f,sh,ci,rn)

# ---- idea 2: sparser classifiers (v2 rule on subset, all circuits measured) ----
import copy
Ktr=d["qkernel"].evaluate(d["X_train"],d["X_sv"])           # train x SV
ftr=Ktr@w+b
def sub_dd(S,w2,b2):
    dd=dict(d); dd["w"]=w2; dd["K_dev"]=K[:,S]; f2=K[:,S]@w2+b2; dd["f_dev"]=f2
    s_star=2.5**2*((w2**2)*K[:,S]*(1-K[:,S])).sum(1)/f2**2
    dd["near"]=s_star>CAP; return dd,f2
def sparse_run(name,params,S,w2,b2,seed=3):
    dd,f2=sub_dd(S,w2,b2)
    f,sh,ci,rn=run(w2,b2,K[:,S],np.arange(len(S)),len(S),1,seed,lazy=False,rounds=L)
    e=evaluate(dd,f,sh)   # strict vs the sparse model's OWN exact classifier
    dis_o=(np.sign(f)!=np.sign(fd)[None]).mean(0)
    r=dict(method=name,params=params,circ=float(ci.mean()),shots=e["mean_shots_per_prediction"],rounds=float(rn.mean()),
      worst=float(e["worst_far"]),strict=e["meets_strict"],worst_vs_orig=float(dis_o[~d["near"]].max()),
      agree_orig=float((np.sign(f)==np.sign(fd)[None]).mean()),acc=float(e["accuracy"]),
      exact_agree_orig=float((np.sign(f2)==np.sign(fd)).mean()),exact_acc=float(((f2>0)==(d["y_dev"]==d["classes"][1])).mean()))
    rows.append(r); print(r,flush=True)
for mm in (16,12,10,8,6):
    S=ordw[:mm]
    sparse_run("top-m |w|, refit b",f"m={mm}",S,w[S],float(np.mean(ftr-Ktr[:,S]@w[S])))
    # reduced-set: LSQ refit of weights+b to original f on training kernel (ridge 1e-3)
    A=np.hstack([Ktr[:,S],np.ones((len(ftr),1))]); sol=np.linalg.solve(A.T@A+1e-3*np.eye(mm+1),A.T@ftr)
    sparse_run("top-m, LSQ refit to f_orig on train",f"m={mm}",S,sol[:-1],sol[-1])
# different C: retrain SVC
y=d["y_train"]; Kfull=d["qkernel"].evaluate(d["X_train"])
for C in (0.1,0.3,0.5,3,10):
    sv2=SVC(kernel="precomputed",C=C).fit(Kfull,y); S=sv2.support_
    # SV set may not be subset of original SVs: evaluate kernel directly
    Kd=d["qkernel"].evaluate(d["X_dev"],d["X_train"][S]); w2=sv2.dual_coef_[0]; b2=sv2.intercept_[0]
    dd=dict(d); dd["w"]=w2; dd["K_dev"]=Kd; f2=Kd@w2+b2; dd["f_dev"]=f2
    dd["near"]=2.5**2*((w2**2)*Kd*(1-Kd)).sum(1)/f2**2>CAP
    f,sh,ci,rn=run(w2,b2,Kd,np.arange(len(S)),len(S),1,4,lazy=False)
    e=evaluate(dd,f,sh); dis_o=(np.sign(f)!=np.sign(fd)[None]).mean(0)
    r=dict(method="retrain C",params=f"C={C},nSV={len(S)}",circ=float(ci.mean()),shots=e["mean_shots_per_prediction"],rounds=float(rn.mean()),
      worst=float(e["worst_far"]),strict=e["meets_strict"],worst_vs_orig=float(dis_o[~d["near"]].max()),
      agree_orig=float((np.sign(f)==np.sign(fd)[None]).mean()),acc=float(e["accuracy"]),
      exact_agree_orig=float((np.sign(f2)==np.sign(fd)).mean()),exact_acc=float(((f2>0)==(d["y_dev"]==d["classes"][1])).mean()))
    rows.append(r); print(r,flush=True)
json.dump(rows,open("rows.json","w"),indent=1)
