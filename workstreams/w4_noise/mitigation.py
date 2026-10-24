"""Reusable mitigation functions for P('00') kernel estimates."""
import numpy as np

def assignment_matrix(cal_counts, shots):
    """cal_counts[prep] = dict bitstring->count for prep in '00','01','10','11' (qiskit order, q1q0).
    Returns A with A[meas, prep] = P(meas|prep)."""
    labels = ['00', '01', '10', '11']
    A = np.zeros((4, 4))
    for j, p in enumerate(labels):
        for i, m in enumerate(labels):
            A[i, j] = cal_counts[p].get(m, 0) / shots
    return A

def ro_mitigate_p00(counts, shots, A):
    """Invert assignment matrix on full 4-outcome distribution; return corrected P('00')
    (unconstrained inverse; clip to [0,1])."""
    labels = ['00', '01', '10', '11']
    p = np.array([counts.get(l, 0) / shots for l in labels])
    q = np.linalg.solve(A, p)
    return float(np.clip(q[0], 0, 1))

def ro_mitigate_p00_vec(P, A):
    """P: (...,4) empirical dists -> corrected P00 array."""
    q = np.linalg.solve(A, P.reshape(-1, 4).T).T
    return np.clip(q[:, 0], 0, 1).reshape(P.shape[:-1])

def fit_affine(k_noisy_cal, k_exact_cal):
    """Least squares k_noisy = a*k_exact + c  -> (a, c)."""
    a, c = np.polyfit(k_exact_cal, k_noisy_cal, 1)
    return a, c

def identity_affine(k_identity_noisy, dim=2):
    """From k(x,x)=1 circuits only, assume full depolarizing toward 1/2^n: k = a*k_true + (1-a)/2^n."""
    k1 = float(np.mean(k_identity_noisy)); floor = 1 / 2**dim
    a = (k1 - floor) / (1 - floor)
    return a, (1 - a) * floor

def correct_kernel(k, a, c):
    return (k - c) / a

def decision(K, w, b, a=1.0, c=0.0):
    """f = sum_i w_i (k_i-c)/a + b  (equivalently f_noisy' = (f_noisy - c*sum w - b)/a + b)."""
    return correct_kernel(K, a, c) @ w + b
