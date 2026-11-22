"""
Phase 6: Global Sensitivity Analysis (Sobol & OAT Sweeps) (Path B)
==================================================================
Executes 14,336 Saltelli model evaluations (N = 1024, D = 6) across parameters:
  - L0 (+/- 10%)
  - K (+/- 25%)
  - alpha (+/- 25%)
  - beta (+/- 25%)
  - eta (0x to 2x)
  - PCSA (+/- 25%)
Using SOBOL_RANDOM_SEED = 42, computing First-Order (S1) and Total-Order (ST)
variance-based sensitivity indices for Table 5 in the manuscript.

Outputs:
  Global_Sensitivity_Analysis_Sobol.xlsx (Table 5) in results/tables/
"""

import os
import glob
import numpy as np
import pandas as pd
import scipy.signal as signal
from scipy import special
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

BASE_DIR = r"d:\Documents\Elbow joint Project"
PATHB_DIR = os.path.join(BASE_DIR, "03_Path_B_Anatomical_Axis_Model")
TABLES_DIR = os.path.join(PATHB_DIR, "results", "tables")
CROSS_DIR = os.path.join(PATHB_DIR, "results", "Cross_Subject_Average")
os.makedirs(TABLES_DIR, exist_ok=True)

OUT_FILE = os.path.join(TABLES_DIR, "Global_Sensitivity_Analysis_Sobol.xlsx")

LIGAMENTS = ['pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'pmCAP', 'plCAP']

LIGAMENT_PROPS = {
    'pRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67, 'L0': 0.02032},
    'LUCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.9953, 'area': 13.35, 'F_uts': 102.4, 'L0': 0.05290},
    'aaUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.81,  'area': 4.31,  'F_uts': 89.2,  'L0': 0.03000},
    'caUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.01,  'area': 4.31,  'F_uts': 89.2,  'L0': 0.02835},
    'paUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 7.35,  'area': 4.31,  'F_uts': 89.2,  'L0': 0.02666},
    'pUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 8.90,  'F_uts': 160.0, 'L0': 0.03019},
    'pmCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.62,  'area': 4.39,  'F_uts': 126.3, 'L0': 0.05310},
    'plCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.96,  'area': 4.39,  'F_uts': 126.3, 'L0': 0.05588}
}

def evaluate_devita_nlsr_vectorized(L_arr, dt, L0_vec, K_vec, alpha_vec, beta_vec, eta_vec, area_vec, F_uts):
    """Vectorized calculation of peak force across thousands of Saltelli parameter sets."""
    # L_arr shape: (T,)
    # Parameter vectors shape: (N,)
    
    # Broadcast strain: (N, T)
    strains = (L_arr[None, :] - L0_vec[:, None]) / L0_vec[:, None]
    
    # Strain rate: (N, T)
    rates = np.gradient(strains, dt, axis=1)
    
    # DeVita computation
    pos = strains > 0
    u = np.zeros_like(strains)
    u[pos] = (strains[pos] / beta_vec[:, None].repeat(strains.shape[1], axis=1)[pos]) ** alpha_vec[:, None].repeat(strains.shape[1], axis=1)[pos]
    
    cdf = np.zeros_like(strains)
    cdf[pos] = 1.0 - np.exp(-u[pos])
    
    shape_p = 1.0 + 1.0 / alpha_vec
    # Gamma functions per parameter set
    first_m = np.zeros_like(strains)
    for i in range(len(alpha_vec)):
        sp = shape_p[i]
        b = beta_vec[i]
        u_row = u[i, :]
        p_row = pos[i, :]
        if np.any(p_row):
            first_m[i, p_row] = b * special.gamma(sp) * special.gammainc(sp, u_row[p_row])
            
    s_elast = np.zeros_like(strains)
    s_elast[pos] = np.maximum(0.0, K_vec[:, None].repeat(strains.shape[1], axis=1)[pos] * (strains[pos] * cdf[pos] - first_m[pos]))
    
    s_visc = np.zeros_like(strains)
    s_visc[pos] = eta_vec[:, None].repeat(strains.shape[1], axis=1)[pos] * rates[pos] * cdf[pos]
    
    s_tot = np.maximum(0.0, s_elast + s_visc)
    f_tot = np.minimum(F_uts, s_tot * area_vec[:, None])
    
    peak_forces = np.max(f_tot, axis=1)
    return peak_forces

def run_gsa():
    print("=" * 70)
    print("RUNNING GLOBAL MONTE CARLO SOBOL SENSITIVITY ANALYSIS (PATH B)")
    print("=" * 70)
    
    # Load Cross-Subject Fast Kinematics
    len_file = os.path.join(CROSS_DIR, "Bicep_CrossSubject_Lengths.xlsx")
    if not os.path.exists(len_file):
        print(f"[ERROR] Cross-subject length file not found: {len_file}")
        return
        
    xls_len = pd.ExcelFile(len_file)
    
    # Problem definition for SALib: 6 parameters (L0, K, alpha, beta, eta, PCSA)
    problem = {
        'num_vars': 6,
        'names': ['L0', 'K', 'alpha', 'beta', 'eta', 'PCSA'],
        'bounds': [
            [0.90, 1.10],  # L0 (+/- 10%)
            [0.75, 1.25],  # K (+/- 25%)
            [0.75, 1.25],  # alpha (+/- 25%)
            [0.75, 1.25],  # beta (+/- 25%)
            [0.00, 2.00],  # eta (0x to 2x)
            [0.75, 1.25]   # PCSA (+/- 25%)
        ]
    }
    
    SOBOL_RANDOM_SEED = 42
    N_SAMPLES = 1024
    param_matrix = sobol_sample.sample(problem, N_SAMPLES, calc_second_order=True, seed=SOBOL_RANDOM_SEED)
    N_TOTAL = param_matrix.shape[0]
    print(f"Generated {N_TOTAL} Saltelli parameter samples across 6 dimensions.")
    
    # Duration for Fast trial: 1.51 s -> 101 points -> dt = 1.51 / 100 = 0.0151 s
    dt = 1.51 / 100.0
    
    sobol_summary_records = []
    sobol_detail_records = []
    
    for lig in LIGAMENTS:
        if lig not in xls_len.sheet_names:
            continue
            
        print(f"\nAnalyzing Global Sobol Indices for: {lig} ...")
        df_l = pd.read_excel(xls_len, lig)
        if 'nor_fast_Mean' not in df_l.columns:
            continue
            
        L_arr = df_l['nor_fast_Mean'].values
        props = LIGAMENT_PROPS[lig]
        
        # Scale parameters from Saltelli matrix
        L0_vec    = param_matrix[:, 0] * props['L0']
        K_vec     = param_matrix[:, 1] * props['K']
        alpha_vec = param_matrix[:, 2] * props['alpha']
        beta_vec  = param_matrix[:, 3] * props['beta']
        eta_vec   = param_matrix[:, 4] * props['eta']
        area_vec  = param_matrix[:, 5] * props['area']
        F_uts     = props['F_uts']
        
        # Evaluate model outputs
        y_forces = evaluate_devita_nlsr_vectorized(L_arr, dt, L0_vec, K_vec, alpha_vec, beta_vec, eta_vec, area_vec, F_uts)
        
        # Sobol Analysis via SALib
        si = sobol_analyze.analyze(problem, y_forces, calc_second_order=True, seed=SOBOL_RANDOM_SEED)
        
        # Record indices
        s1_l0, st_l0 = float(si['S1'][0]), float(si['ST'][0])
        s1_k,  st_k  = float(si['S1'][1]), float(si['ST'][1])
        s1_a,  st_a  = float(si['S1'][2]), float(si['ST'][2])
        s1_b,  st_b  = float(si['S1'][3]), float(si['ST'][3])
        s1_e,  st_e  = float(si['S1'][4]), float(si['ST'][4])
        s1_p,  st_p  = float(si['S1'][5]), float(si['ST'][5])
        
        # Interactions = ST - S1
        int_l0 = max(0.0, st_l0 - s1_l0)
        
        sobol_summary_records.append({
            'Ligament': lig,
            'Baseline_Force_N': round(float(np.mean(y_forces)), 2),
            'S1_L0': round(s1_l0, 4),
            'ST_L0': round(st_l0, 4),
            'S1_K': round(s1_k, 4),
            'ST_K': round(st_k, 4),
            'S1_alpha': round(s1_a, 4),
            'ST_alpha': round(st_a, 4),
            'S1_beta': round(s1_b, 4),
            'ST_beta': round(st_b, 4),
            'S1_eta': round(s1_e, 4),
            'ST_eta': round(st_e, 4),
            'S1_PCSA': round(s1_p, 4),
            'ST_PCSA': round(st_p, 4),
            'L0_Total_Variance_Pct': f"{st_l0 * 100.0:.1f}%",
            'Interaction_Index_L0': round(int_l0, 4)
        })
        
        for p_idx, p_name in enumerate(problem['names']):
            sobol_detail_records.append({
                'Ligament': lig,
                'Parameter': p_name,
                'S1': float(si['S1'][p_idx]),
                'S1_conf': float(si['S1_conf'][p_idx]),
                'ST': float(si['ST'][p_idx]),
                'ST_conf': float(si['ST_conf'][p_idx])
            })
            
    df_sobol_sum = pd.DataFrame(sobol_summary_records)
    df_sobol_det = pd.DataFrame(sobol_detail_records)
    
    with pd.ExcelWriter(OUT_FILE) as writer:
        df_sobol_sum.to_excel(writer, sheet_name="Table5_Sobol_Summary", index=False)
        df_sobol_det.to_excel(writer, sheet_name="Full_Sobol_Indices", index=False)
        
    print(f"\n[SUCCESS] Table 5 Global Sensitivity Analysis saved to: {OUT_FILE}")
    print(df_sobol_sum[['Ligament', 'S1_L0', 'ST_L0', 'S1_K', 'ST_K', 'S1_eta', 'ST_eta', 'L0_Total_Variance_Pct']])

if __name__ == "__main__":
    run_gsa()
