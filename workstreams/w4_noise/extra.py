import sys, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness'))
from qhack_harness import load
d=load(); w,b,f=d['w'],d['b'],d['f_dev']
print('sum w = %.3e, b=%.4f, n_dev=%d'%(w.sum(),b,len(f)))
# model from measured fit (analyze_out): raw a=0.8877,c=0.0282 -> f' = a(f-b)+b + c*sum(w)
for a in (0.8877,0.9161):
    fn=a*(f-b)+b
    print('a=%.4f model-based (full dev): disagreements %d / %d; threshold shift b(1-a)=%.4f; near=%d, disagreers near=%s'%(a,(np.sign(fn)!=np.sign(f)).sum(),len(f),b*(1-a),d['near'].sum(),d['near'][np.sign(fn)!=np.sign(f)].tolist()))
