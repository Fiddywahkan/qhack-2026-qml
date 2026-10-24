from sim import *
d = load(); print('N dev',len(d['f_dev']),'near',d['near'].sum())
def show(r): print(f"{r['rule']:10s} {str(r['params']):28s} {r['kind']:9s} mean={r['mean']:.0f} p95={r['p95']:.0f} rnd={r['rounds']:.2f} worst={r['worst']:.3f} strict={r['strict']} agree={r['agree']:.4f} acc={r['acc']:.4f} cappts={r['cap_pts']}")
for dl in (0.05,0.1,0.2,0.4,0.8,1.5):
    show(row('bonf',f'delta={dl}','tuned' if dl>0.05 else 'principled',d,zbonf(dl)))
for rho in (1,10,100,1000):
    show(row('mix',f'd=.05 rho={rho}','principled',d,zmix(0.05,rho)))
show(row('poly','d=.05 p=2','principled',d,zpoly(0.05)))
