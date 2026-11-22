import os, glob
import numpy as np
import pandas as pd
from scipy import special, signal, optimize

cross_len_path = r"PROJECT_DIR\results\Cross_Subject_Average\Bicep_CrossSubject_Lengths.xlsx"
xf = pd.ExcelFile(cross_len_path)

dt = 0.01

def compute_force_series(L, dt, L0, alpha, beta, K, area, eta_val, F_uts=110.0, fc=8.0):
    wn = fc / (0.5 * (1.0/dt))
    b, a = signal.butter(4, wn, btype='low')
    L_f = signal.filtfilt(b, a, L)
    v_f = np.gradient(L_f, dt)
    
    eps = (L_f - L0) / L0
    epsdot = v_f / L0
    
    pos = eps > 0
    u = np.zeros_like(eps)
    u[pos] = (eps[pos] / beta)**alpha
    cdf = np.zeros_like(eps)
    cdf[pos] = 1.0 - np.exp(-u[pos])
    shape_p = 1.0 + 1.0 / alpha
    first_m = np.zeros_like(eps)
    first_m[pos] = beta * special.gamma(shape_p) * special.gammainc(shape_p, u[pos])
    
    s_elast = np.zeros_like(eps)
    s_elast[pos] = np.maximum(0.0, K * (eps[pos] * cdf[pos] - first_m[pos]))
    s_visc = np.zeros_like(eps)
    s_visc[pos] = eta_val * epsdot[pos] * cdf[pos]
    
    s_tot = np.maximum(0.0, s_elast + s_visc)
    F_tot = np.minimum(F_uts, s_tot * area)
    F_elast = np.minimum(F_uts, s_elast * area)
    F_visc = s_visc * area
    
    return {'F_tot': F_tot, 'F_elast': F_elast, 'F_visc': F_visc, 'eps': eps, 'epsdot': epsdot}

df = pd.read_excel(xf, 'pRCL')
cols = [c for c in df.columns if c != 'Cycle %']
print("Columns in pRCL:", cols)
fast_cols = [c for c in cols if 'fast' in c.lower()]
L = df[fast_cols].mean(axis=1).values if fast_cols else df[cols].mean(axis=1).values

L0 = 0.02032
K = 49.2
area = 5.67
alpha = 1.01
beta = 0.080

for eta in [0.01, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0]:
    res = compute_force_series(L, dt, L0, alpha, beta, K, area, eta)
    idx_peak = np.argmax(res['F_tot'])
    ft = res['F_tot'][idx_peak]
    fv = res['F_visc'][idx_peak]
    fe = res['F_elast'][idx_peak]
    sh = fv / ft if ft > 0 else 0
    print(f"eta={eta:5.2f}: F_tot={ft:5.2f} N, F_elast={fe:5.2f} N, F_visc={fv:5.2f} N, share={sh*100:5.1f}%, eps={res['eps'][idx_peak]*100:.2f}%, epsdot={res['epsdot'][idx_peak]:.3f}")
