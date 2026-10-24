import sys,numpy as np
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'harness'))
from qhack_harness import load
K=load(2)["K_dev"]; print("shot-noise floor mean|dk| at S=4096: %.4f"%(0.798*np.sqrt(K*(1-K)/4096)).mean())
