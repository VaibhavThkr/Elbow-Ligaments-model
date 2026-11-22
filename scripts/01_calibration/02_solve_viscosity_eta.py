"""
Recalibrate eta* for aRCL, pRCL, and LUCL using Brent's method
to enforce exactly 30% viscous share at peak force in Fast tempo (1.0 s cycle).
"""
import os, glob
import numpy as np
import pandas as pd
from scipy import special, signal, optimize

# Fast tempo length trajectory on cohort mean
# Let's load the cohort average length from Bicep_CrossSubject_Lengths.xlsx
cross_len_path = r"PROJECT_DIR\results\Cross_Subject_Average\Bicep_CrossSubject_Lengths.xlsx"
xf = pd.ExcelFile(cross_len_path)

dt = 0.01 # 100 Hz time-normalized or 1.0s total duration for Fast tempo
# In Fast tempo, movement duration is 1.0s (101 points, dt = 0.01s)

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

def solve_eta_for_ligament(lig_name, L0, alpha, beta, K, area, F_uts, target_share=0.30):
    df = pd.read_excel(xf, lig_name)
    # fast column: usually '1B_nor_fast_right' or 'nor_fast' or fast column
    fast_cols = [c for c in df.columns if 'fast' in c.lower()]
    if fast_cols:
        L = df[fast_cols].mean(axis=1).values
    else:
        cols = [c for c in df.columns if c != 'Cycle %']
        L = df[cols].mean(axis=1).values
        
    def objective(eta_cand):
        res = compute_force_series(L, dt, L0, alpha, beta, K, area, eta_cand, F_uts)
        idx_peak = np.argmax(res['F_tot'])
        if res['F_tot'][idx_peak] <= 1e-4:
            return -target_share
        actual_share = res['F_visc'][idx_peak] / res['F_tot'][idx_peak]
        return actual_share - target_share
        
    eta_opt = optimize.brentq(objective, 0.01, 50.0)
    res_opt = compute_force_series(L, dt, L0, alpha, beta, K, area, eta_opt, F_uts)
    idx_peak = np.argmax(res_opt['F_tot'])
    
    print(f"=== {lig_name} ===")
    print(f"  Calibrated eta*: {eta_opt:.4f} MPa*s")
    print(f"  Peak Total Force: {res_opt['F_tot'][idx_peak]:.2f} N")
    print(f"  Peak Elastic Force: {res_opt['F_elast'][idx_peak]:.2f} N")
    print(f"  Peak Viscous Force: {res_opt['F_visc'][idx_peak]:.2f} N")
    print(f"  Viscous Share at Peak: {res_opt['F_visc'][idx_peak]/res_opt['F_tot'][idx_peak]*100:.1f}%")
    print(f"  Peak Strain: {res_opt['eps'][idx_peak]*100:.2f}%")
    print(f"  Peak Strain Rate: {res_opt['epsdot'][idx_peak]:.3f} s^-1")
    return eta_opt

print("Solving eta* on Cohort Mean Fast Trajectory:")
# pRCL: area = 5.67 mm^2, K = 49.2 MPa, alpha=1.01, beta=0.080, L0 ~ 0.02032 m
eta_pRCL = solve_eta_for_ligament('pRCL', 0.02032, 1.01, 0.080, 49.2, 5.67, 36.67)

# LUCL: area = 13.35 mm^2, K = 49.2 MPa, alpha=1.01, beta=0.080, L0 ~ 0.05290 m
eta_LUCL = solve_eta_for_ligament('LUCL', 0.05290, 1.01, 0.080, 49.2, 13.35, 102.4)

# aRCL: area = 5.67 mm^2, K = 49.2 MPa, alpha=1.01, beta=0.080, L0 ~ 0.02033 m
eta_aRCL = solve_eta_for_ligament('aRCL', 0.02033, 1.01, 0.080, 49.2, 5.67, 36.67)
