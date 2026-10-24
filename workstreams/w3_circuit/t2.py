import sys, numpy as np
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness'))
from qhack_harness import load
print("reps | dev acc (sign f vs y) | #SV | oracle mean min(s*,2048) | #near-tie | agree w/ reps2? n/a")
for r in (1,2,3):
    d=load(r); K,w,f=d["K_dev"],d["w"],d["f_dev"]
    s=2.5**2*((w**2)*K*(1-K)).sum(1)/f**2
    pred=np.where(f>0,d["classes"][1],d["classes"][0])
    acc=(pred==d["y_dev"]).mean()
    print(r,"acc %.4f nSV %d shots/circ(oracle) %.1f  per-pred %.0f near %d  ndev %d"%(acc,len(w),np.minimum(s,2048).mean(),np.minimum(s,2048).mean()*len(w),d["near"].sum(),len(f)))
    print("   train acc %.3f"%d["svm"].score(d["qkernel"].evaluate(d["X_train"]),d["y_train"]))
