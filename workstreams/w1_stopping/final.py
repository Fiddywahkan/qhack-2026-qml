import sys, json
from sim import *
import calib
d = load(); T=int(sys.argv[1]) if len(sys.argv)>1 else 200
SEED=0
obf_sp=lambda t,a:2*(1-norm.cdf(norm.isf(a/2)/np.sqrt(t)))
rules=[]
def add(name,params,kind,zf,ney=False): rules.append((name,params,kind,zf,ney))
add('v2 Bonferroni','delta=.05,L=9','principled',zbonf(.05))
add('Bonferroni','delta=.1','tuned',zbonf(.1))
add('Bonferroni','delta=.2','tuned',zbonf(.2))
add('Bonferroni','delta=.3','tuned',zbonf(.3))
for dl in (.05,.1,.2):
    zp=calib.pocock(dl); add('Pocock GS (MC-calibrated const z)',f'delta={dl},z={zp:.3f}','principled' if dl==.05 else 'tuned',zspend([norm.sf(zp)]*L) if False else (lambda zp=zp:(lambda j,s:zp))())
add('Lan-DeMets Pocock-type spend','delta=.05','principled',zspend(norm.sf(calib.spend_bounds(.05,calib.pock1))))
add('Lan-DeMets Pocock-type spend','delta=.1','tuned',zspend(norm.sf(calib.spend_bounds(.1,calib.pock1))))
add('Lan-DeMets OBF-type spend','delta=.05','principled',zspend(norm.sf(calib.spend_bounds(.05,obf_sp))))
add('Lan-DeMets OBF-type spend','delta=.1','tuned',zspend(norm.sf(calib.spend_bounds(.1,obf_sp))))
add('poly spend p=2','delta=.05','principled',zpoly(.05))
add('Robbins normal-mixture CS','delta=.05,rho=1','principled',zmix(.05,1))
base=[r for r in rules]
for r in base:
    if r[0] in('v2 Bonferroni','Bonferroni','Pocock GS (MC-calibrated const z)','Lan-DeMets Pocock-type spend'):
        add(r[0]+' + Neyman',r[1],r[2],r[3],True)
out=[]
hdr="| rule | params | kind | mean shots/pred | p95 | mean rounds | worst non-near disagr. | strict | agree | dev acc | cap pts(>50%) | rounds dist 1..9 |"
print(hdr)
for n,p,k,zf,ney in rules:
    r=row(n,p,k,d,zf,ney,T,SEED); out.append({**r,'rounds_hist':r['rounds_hist'].tolist()})
    print(f"| {n} | {p} | {k} | {r['mean']:.0f} | {r['p95']:.0f} | {r['rounds']:.2f} | {r['worst']:.3f} | {r['strict']} | {r['agree']:.4f} | {r['acc']:.4f} | {r['cap_pts']}/{r['near']} near | {np.round(r['rounds_hist'],2).tolist()} |")
json.dump(out,open(f'results_T{T}.json','w'),default=float)
