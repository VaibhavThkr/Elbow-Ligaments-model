import pandas as pd
import numpy as np
from scipy import special, signal, optimize

cross_len_path = r"PROJECT_DIR\results\Cross_Subject_Average\Bicep_CrossSubject_Lengths.xlsx"
xf = pd.ExcelFile(cross_len_path)
dt = 0.01

def compute_force(L, dt, L0, alpha, beta, K, area, eta, fc=8.0):
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
    s_visc[pos] = eta * epsdot[pos] * cdf[pos]
    
    s_tot = np.maximum(0.0, s_elast + s_visc)
    F_tot = s_tot * area
    F_elast = s_elast * area
    F_visc = s_visc * area
    return F_tot, F_elast, F_visc, eps, epsdot

def solve_eta(lig_name, L0, area, K=49.2, alpha=1.01, beta=0.080):
    df = pd.read_excel(xf, lig_name)
    L = df['nor_fast_Mean'].values
    
    def objective(eta_cand):
        Ft, Fe, Fv, eps, edot = compute_force(L, dt, L0, alpha, beta, K, area, eta_cand)
        idx = np.argmax(Ft)
        return (Fv[idx] / Ft[idx]) - 0.30
        
    eta_opt = optimize.brentq(objective, 0.1, 20.0)
    Ft, Fe, Fv, eps, edot = compute_force(L, dt, L0, alpha, beta, K, area, eta_opt)
    idx = np.argmax(Ft)
    print(f"=== {lig_name} ===")
    print(f"  Optimized eta*: {eta_opt:.4f} MPa*s")
    print(f"  Peak F_tot: {Ft[idx]:.2f} N (Elastic: {Fe[idx]:.2f} N, Viscous: {Fv[idx]:.2f} N)")
    print(f"  Viscous share: {Fv[idx]/Ft[idx]*100:.2f}%")
    print(f"  Peak strain: {eps[idx]*100:.2f}%")
    print(f"  Peak strain rate: {edot[idx]:.3f} s^-1")
    return eta_opt

print("=== Solving exact 30% viscous share ===")
eta_pRCL = solve_eta('pRCL', 0.02030, 5.67)
eta_LUCL = solve_eta('LUCL', 0.05290, 13.35)
