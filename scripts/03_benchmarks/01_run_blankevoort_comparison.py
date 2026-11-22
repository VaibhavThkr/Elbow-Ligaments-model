"""
Phase 5 (Part 2): Method A vs Method B Benchmark Comparison (Table 6) (Anatomical Axis Model)
==============================================================================
Compares Elastic Blankevoort (1991) vs. Damped Blankevoort (OpenSim 4.x) vs.
Calibrated DeVita NLSR (Anatomical Axis Elbow Model) across all 6 participants.

Outputs:
  Table6_MethodA_Model_Comparison_Benchmark.xlsx in results/tables/
"""

import os
import glob
import re
import numpy as np
import pandas as pd
import scipy.special as special

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
MODEL_DIR = os.path.join(PROJECT_DIR, "models", "calibrated_devita_nlsr")
SIM_DIR = os.path.join(PROJECT_DIR, "simulations")
TABLES_DIR = os.path.join(PROJECT_DIR, "results", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

OUT_FILE = os.path.join(TABLES_DIR, "Table6_MethodA_Model_Comparison_Benchmark.xlsx")

SUBJECTS = [
    ("S1", "Gokul",  "1_S1 Bicep"),
    ("S2", "Sanyam", "2_S2 Bicep"),
    ("S3", "Sunil",  "3_S3 Bicep"),
    ("S4", "Piyush", "4_S4 Bicep"),
    ("S5", "Nitin",  "5_S5 Bicep"),
    ("S6", "Rahul",  "6_S6 Bicep")
]

PRIMARY_LIGAMENTS = ['pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'pmCAP', 'plCAP']

LIGAMENT_PROPS = {
    'pRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67, 'eps_t': 0.03, 'k_lin': 85.46,  'c_damp': 9.87},
    'LUCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.9953, 'area': 13.35, 'F_uts': 102.4, 'eps_t': 0.03, 'k_lin': 212.85, 'c_damp': 22.60},
    'aaUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.81,  'area': 4.31,  'F_uts': 89.2,  'eps_t': 0.03, 'c_damp': 15.0},
    'caUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.01,  'area': 4.31,  'F_uts': 89.2,  'eps_t': 0.03, 'c_damp': 15.0},
    'paUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 7.35,  'area': 4.31,  'F_uts': 89.2,  'eps_t': 0.03, 'c_damp': 15.0},
    'pUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 8.90,  'F_uts': 160.0, 'eps_t': 0.03, 'c_damp': 25.0},
    'pmCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.62,  'area': 4.39,  'F_uts': 126.3, 'eps_t': 0.03, 'c_damp': 20.0},
    'plCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.96,  'area': 4.39,  'F_uts': 126.3, 'eps_t': 0.03, 'c_damp': 20.0}
}

TEMPOS = [('nor_fast', 'Fast'), ('nor_med', 'Medium'), ('nor_slow', 'Slow')]

def f_blankevoort_elastic(eps, k_lin, eps_t):
    f = np.zeros_like(eps)
    m1 = (eps > 0) & (eps <= 2.0 * eps_t)
    m2 = eps > 2.0 * eps_t
    f[m1] = 0.5 * k_lin * (eps[m1]**2) / eps_t
    f[m2] = k_lin * (eps[m2] - eps_t)
    return np.maximum(0.0, f)

def f_blankevoort_damped(eps, epsdot, k_lin, eps_t, c_lin):
    f_spring = f_blankevoort_elastic(eps, k_lin, eps_t)
    mask = eps > 0
    f_damp = np.zeros_like(eps)
    f_damp[mask] = c_lin * epsdot[mask]
    f_tot = np.zeros_like(eps)
    f_tot[mask] = np.maximum(0.0, f_spring[mask] + f_damp[mask])
    return f_tot, f_damp

def f_devita_nlsr(eps, epsdot, props):
    alpha = props['alpha']
    beta  = props['beta']
    K     = props['K']
    eta   = props['eta']
    area  = props['area']
    F_uts = props['F_uts']
    
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
    F_tot = np.minimum(F_uts, s_tot * area)
    F_visc = np.minimum(F_uts, s_visc * area)
    F_elast = np.minimum(F_uts, s_elast * area)
    return F_tot, F_visc, F_elast

def get_calibrated_l0_from_osim(model_path):
    with open(model_path, 'r', encoding='utf-8') as f:
        text = f.read()
    l0_dict = {}
    for lig in PRIMARY_LIGAMENTS:
        m = re.search(r'<NonLinearLigament name="' + lig + r'">.*?</NonLinearLigament>', text, re.DOTALL)
        if m:
            block = m.group(0)
            rl_match = re.search(r'<resting_length>(.*?)</resting_length>', block)
            if rl_match:
                l0_dict[lig] = float(rl_match.group(1))
    if 'RCL' not in l0_dict: l0_dict['RCL'] = 0.01981
    if 'LUCL' not in l0_dict: l0_dict['LUCL'] = 0.05250
    return l0_dict

def run_comparison():
    print("=" * 70)
    print("RUNNING TABLE 6: METHOD A VS METHOD B BENCHMARK COMPARISON (ANATOMICAL AXIS MODEL)")
    print("=" * 70)
    
    # Preload length and strain rate workbooks for all subjects
    subject_data = {}
    for sid, name, folder_name in SUBJECTS:
        ar_dir = os.path.join(SIM_DIR, folder_name, "Analysis Results")
        f_len = os.path.join(ar_dir, f"{name} Bicep_Lengths.xlsx")
        f_rates = os.path.join(ar_dir, f"{name} Bicep_StrainRates.xlsx")
        f_strains = os.path.join(ar_dir, f"{name} Bicep_Strains.xlsx")
        model_path = os.path.join(MODEL_DIR, f"Elbow_model-{sid}.osim")
        
        if os.path.exists(f_len) and os.path.exists(f_rates):
            subject_data[sid] = {
                'lengths': pd.ExcelFile(f_len),
                'rates': pd.ExcelFile(f_rates),
                'strains': pd.ExcelFile(f_strains),
                'l0_map': get_calibrated_l0_from_osim(model_path)
            }
            
    records = []
    
    for lig in PRIMARY_LIGAMENTS:
        props = LIGAMENT_PROPS[lig]
        k_lin = props.get('k_lin', props['K'] * props['area'])
        eps_t = props['eps_t']
        c_damp = props['c_damp']
        
        for cond_col, tempo_name in TEMPOS:
            devita_peaks = []
            elastic_peaks = []
            damped_peaks = []
            devita_visc_shares = []
            damped_visc_shares = []
            
            for sid, sdata in subject_data.items():
                if lig in sdata['strains'].sheet_names:
                    df_str = pd.read_excel(sdata['strains'], lig)
                    df_rate = pd.read_excel(sdata['rates'], lig)
                    
                    if cond_col in df_str.columns and cond_col in df_rate.columns:
                        eps = df_str[cond_col].values
                        epsdot = df_rate[cond_col].values
                        
                        # DeVita NLSR
                        f_dev, f_dev_v, f_dev_e = f_devita_nlsr(eps, epsdot, props)
                        pk_dev = np.max(f_dev)
                        idx_dev = np.argmax(f_dev)
                        v_share_dev = (f_dev_v[idx_dev] / pk_dev * 100.0) if pk_dev > 0.01 else 0.0
                        
                        # Blankevoort Elastic
                        f_elast = f_blankevoort_elastic(eps, k_lin, eps_t)
                        pk_elast = np.max(f_elast)
                        
                        # Blankevoort Damped
                        f_damp_tot, f_damp_v = f_blankevoort_damped(eps, epsdot, k_lin, eps_t, c_damp)
                        pk_damped = np.max(f_damp_tot)
                        idx_damp = np.argmax(f_damp_tot)
                        v_share_damp = (f_damp_v[idx_damp] / pk_damped * 100.0) if pk_damped > 0.01 else 0.0
                        
                        devita_peaks.append(pk_dev)
                        elastic_peaks.append(pk_elast)
                        damped_peaks.append(pk_damped)
                        devita_visc_shares.append(v_share_dev)
                        damped_visc_shares.append(v_share_damp)
                        
            if devita_peaks:
                m_dev, s_dev = np.mean(devita_peaks), np.std(devita_peaks, ddof=1)
                m_ela, s_ela = np.mean(elastic_peaks), np.std(elastic_peaks, ddof=1)
                m_dmp, s_dmp = np.mean(damped_peaks), np.std(damped_peaks, ddof=1)
                m_v_dev = np.mean(devita_visc_shares)
                m_v_dmp = np.mean(damped_visc_shares)
                
                # Difference Method B vs Method A Elastic
                diff_pct_ela = ((m_dev - m_ela) / m_ela * 100.0) if m_ela > 0.01 else 0.0
                diff_pct_dmp = ((m_dev - m_dmp) / m_dmp * 100.0) if m_dmp > 0.01 else 0.0
                
                records.append({
                    'Ligament': lig,
                    'Tempo': tempo_name,
                    'MethodA_Elastic_Peak_N': f"{m_ela:.2f} ± {s_ela:.2f}",
                    'MethodA_Damped_Peak_N': f"{m_dmp:.2f} ± {s_dmp:.2f}",
                    'MethodB_DeVita_NLSR_Peak_N': f"{m_dev:.2f} ± {s_dev:.2f}",
                    'MethodA_Damped_Viscous_Share_Pct': f"{m_v_dmp:.1f}%",
                    'MethodB_DeVita_Viscous_Share_Pct': f"{m_v_dev:.1f}%",
                    'Diff_DeVita_vs_Elastic_Pct': f"{diff_pct_ela:+.1f}%",
                    'Diff_DeVita_vs_Damped_Pct': f"{diff_pct_dmp:+.1f}%",
                    'Raw_DeVita_Mean': m_dev,
                    'Raw_Elastic_Mean': m_ela,
                    'Raw_Damped_Mean': m_dmp
                })
                
    df_comp = pd.DataFrame(records)
    with pd.ExcelWriter(OUT_FILE) as writer:
        df_comp.to_excel(writer, sheet_name="Table6_Model_Comparison", index=False)
        
    print(f"\n[SUCCESS] Table 6 benchmark report saved to: {OUT_FILE}")
    print(df_comp[['Ligament', 'Tempo', 'MethodA_Elastic_Peak_N', 'MethodA_Damped_Peak_N', 'MethodB_DeVita_NLSR_Peak_N', 'Diff_DeVita_vs_Elastic_Pct']].head(12))

if __name__ == "__main__":
    run_comparison()
