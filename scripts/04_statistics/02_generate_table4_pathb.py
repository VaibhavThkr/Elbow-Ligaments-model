"""
Phase 4 (Part 2): Generate Table 4 & Participant Peak Metrics (Path B)
=====================================================================
Processes the generated participant Excel workbooks from Path B to produce:
  1. Table4_Participant_and_Cohort_Peaks.xlsx in results/tables/
  2. Calibrated_Model_Participant_Level_Results_and_Repetitions.xlsx in results/tables/

Integrates Option A: 3-Fiber RCL Fan (aRCL, cRCL, pRCL, RCL_Total) + LUCL + Medial + Capsular.
"""

import os
import glob
import re
import numpy as np
import pandas as pd
import scipy.special as special
from scipy.signal import find_peaks

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

BASE_DIR = r"d:\Documents\Elbow joint Project"
PATHB_DIR = os.path.join(BASE_DIR, "03_Path_B_Anatomical_Axis_Model")
MODEL_DIR = os.path.join(PATHB_DIR, "models", "calibrated_devita_nlsr")
SIM_DIR = os.path.join(PATHB_DIR, "simulations")
TABLES_DIR = os.path.join(PATHB_DIR, "results", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

OUT_FILE = os.path.join(TABLES_DIR, "Table4_Participant_and_Cohort_Peaks.xlsx")
COMPREHENSIVE_OUT = os.path.join(TABLES_DIR, "Calibrated_Model_Participant_Level_Results_and_Repetitions.xlsx")

SUBJECTS = [
    ("S1", "Gokul",  "1_S1 Bicep", "1B"),
    ("S2", "Sanyam", "2_S2 Bicep", "2B"),
    ("S3", "Sunil",  "3_S3 Bicep", "3B"),
    ("S4", "Piyush", "4_S4 Bicep", "4B"),
    ("S5", "Nitin",  "5_S5 Bicep", "5B"),
    ("S6", "Rahul",  "6_S6 Bicep", "6B")
]

PRIMARY_LIGAMENTS = ['pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'pmCAP', 'plCAP']
ALL_LIGAMENTS = [
    'pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'tUCL',
    'pmCAP', 'plCAP', 'amCAP', 'alCAP',
    'aAL', 'pAL', 'ALCL', 'aRCL', 'cRCL', 'RCL', 'QL', 'OC'
]

ANATOMICAL_GROUP = {
    'pRCL': 'Lateral Collateral (Posterior Functional Bundle)',
    'cRCL': 'Lateral Collateral (Central Isometric Axis Bundle)',
    'aRCL': 'Lateral Collateral (Anterior Bundle)',
    'RCL': 'Lateral Collateral (Composite Fan Total)',
    'LUCL': 'Lateral Ulnar Collateral',
    'aAL': 'Anterior Annular',
    'pAL': 'Posterior Annular',
    'ALCL': 'Lateral Accessory',
    'aaUCL': 'Medial Collateral (Anterior-Anterior)',
    'caUCL': 'Medial Collateral (Central-Anterior)',
    'paUCL': 'Medial Collateral (Posterior-Anterior)',
    'pUCL': 'Medial Collateral (Posterior Bundle)',
    'tUCL': 'Medial Collateral (Transverse)',
    'amCAP': 'Anterior Capsule (Medial)',
    'alCAP': 'Anterior Capsule (Lateral)',
    'pmCAP': 'Posterior Capsule (Medial)',
    'plCAP': 'Posterior Capsule (Lateral)',
    'QL': 'Quadrate Ligament',
    'OC': 'Oblique Cord'
}

LIGAMENT_PROPS = {
    'aRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'cRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'pRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'RCL':   {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 17.01, 'F_uts': 110.0},
    'LUCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.9953, 'area': 13.35, 'F_uts': 102.4},
    'aaUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.81,  'area': 4.31,  'F_uts': 89.2},
    'caUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.01,  'area': 4.31,  'F_uts': 89.2},
    'paUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 7.35,  'area': 4.31,  'F_uts': 89.2},
    'pUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 8.90,  'F_uts': 160.0},
    'tUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 2.23,  'F_uts': 160.0},
    'pmCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.62,  'area': 4.39,  'F_uts': 126.3},
    'plCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.96,  'area': 4.39,  'F_uts': 126.3},
    'amCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.62,  'area': 4.39,  'F_uts': 126.3},
    'alCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.96,  'area': 4.39,  'F_uts': 126.3},
    'aAL':   {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.94,  'area': 6.48,  'F_uts': 102.95},
    'pAL':   {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 9.87,  'area': 6.48,  'F_uts': 102.95},
    'ALCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 2.44,  'area': 4.25,  'F_uts': 65.0},
    'QL':    {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 0.00,  'area': 0.00,  'F_uts': 0.0},
    'OC':    {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 0.00,  'area': 0.00,  'F_uts': 0.0}
}

CONDITION_MAP = {
    'nor_fast': ('Normal', 'Fast', 'Normal Fast'),
    'nor_med':  ('Normal', 'Medium', 'Normal Medium'),
    'nor_slow': ('Normal', 'Slow', 'Normal Slow'),
    'rest_fast': ('Restricted', 'Fast', 'Restricted Fast'),
    'rest_med':  ('Restricted', 'Medium', 'Restricted Medium'),
    'rest_slow': ('Restricted', 'Slow', 'Restricted Slow')
}

def get_calibrated_l0_from_osim(model_path):
    with open(model_path, 'r', encoding='utf-8') as f:
        text = f.read()
    l0_dict = {}
    for lig in ALL_LIGAMENTS:
        m = re.search(r'<NonLinearLigament name="' + lig + r'">.*?</NonLinearLigament>', text, re.DOTALL)
        if m:
            block = m.group(0)
            rl_match = re.search(r'<resting_length>(.*?)</resting_length>', block)
            if rl_match:
                l0_dict[lig] = float(rl_match.group(1))
    if 'aRCL' not in l0_dict: l0_dict['aRCL'] = 0.02000
    if 'cRCL' not in l0_dict: l0_dict['cRCL'] = 0.01973
    if 'pRCL' not in l0_dict: l0_dict['pRCL'] = 0.02030
    if 'LUCL' not in l0_dict: l0_dict['LUCL'] = 0.05250
    if 'RCL'  not in l0_dict: l0_dict['RCL']  = 0.01981
    return l0_dict

def extract_participant_peaks():
    print("=" * 70)
    print("EXTRACTING PATH B PARTICIPANT-LEVEL PEAKS & LIGAMENT METRICS")
    print("=" * 70)
    
    rows = []
    for sid, name, folder_name, prefix in SUBJECTS:
        ar_dir = os.path.join(SIM_DIR, folder_name, "Analysis Results")
        f_forces  = os.path.join(ar_dir, f"{name} Bicep_Forces.xlsx")
        f_strains = os.path.join(ar_dir, f"{name} Bicep_Strains.xlsx")
        f_rates   = os.path.join(ar_dir, f"{name} Bicep_StrainRates.xlsx")
        f_lengths = os.path.join(ar_dir, f"{name} Bicep_Lengths.xlsx")
        
        if not (os.path.exists(f_forces) and os.path.exists(f_strains) and os.path.exists(f_rates)):
            print(f"[SKIP] Missing Excel workbooks for {sid} ({name})")
            continue
            
        model_path = os.path.join(MODEL_DIR, f"Elbow_model-{sid}.osim")
        l0_map = get_calibrated_l0_from_osim(model_path)
        
        xls_f = pd.ExcelFile(f_forces)
        xls_s = pd.ExcelFile(f_strains)
        xls_r = pd.ExcelFile(f_rates)
        xls_l = pd.ExcelFile(f_lengths)
        
        for lig in ALL_LIGAMENTS:
            if lig not in xls_f.sheet_names:
                continue
            props = LIGAMENT_PROPS[lig]
            k_lin = props['K'] * props['area']
            l0_val = l0_map.get(lig, 0.02) * 1000.0
            
            df_f = pd.read_excel(xls_f, lig)
            df_s = pd.read_excel(xls_s, lig)
            df_r = pd.read_excel(xls_r, lig)
            df_l = pd.read_excel(xls_l, lig)
            
            cyc = df_f['Cycle %'].values
            cond_cols = [c for c in df_f.columns if c != 'Cycle %']
            
            for col in cond_cols:
                forces  = df_f[col].values
                strains = df_s[col].values
                rates   = df_r[col].values
                lengths = df_l[col].values * 1000.0
                
                cond_info = CONDITION_MAP.get(col, ('Unknown', 'Unknown', col))
                cond_type, tempo, cond_label = cond_info
                
                pk_force = np.max(forces)
                mn_force = np.mean(forces)
                pk_idx   = np.argmax(forces)
                timing_pk = cyc[pk_idx]
                
                pk_strain = np.max(strains) * 100.0
                mn_strain = np.mean(strains) * 100.0
                pk_rate   = np.max(rates)
                
                # DeVita elastic and viscous force components
                alpha, beta, K, eta, area, F_uts = props['alpha'], props['beta'], props['K'], props['eta'], props['area'], props['F_uts']
                pos = strains > 0
                u = np.zeros_like(strains)
                u[pos] = (strains[pos] / beta)**alpha
                cdf = np.zeros_like(strains)
                cdf[pos] = 1.0 - np.exp(-u[pos])
                shape_p = 1.0 + 1.0 / alpha
                first_m = np.zeros_like(strains)
                first_m[pos] = beta * special.gamma(shape_p) * special.gammainc(shape_p, u[pos])
                s_elast = np.zeros_like(strains)
                s_elast[pos] = np.maximum(0.0, K * (strains[pos] * cdf[pos] - first_m[pos]))
                s_visc = np.zeros_like(strains)
                s_visc[pos] = eta * rates[pos] * cdf[pos]
                
                f_elast = np.minimum(F_uts, s_elast * area)
                f_visc  = np.minimum(F_uts, s_visc * area)
                
                pk_elast = np.max(f_elast)
                pk_visc  = np.max(f_visc)
                visc_at_pk = f_visc[pk_idx]
                visc_share = (visc_at_pk / pk_force) * 100.0 if pk_force > 0.01 else 0.0
                
                rows.append({
                    'Participant': sid,
                    'Participant_Name': name,
                    'Condition': col,
                    'Condition_Type': cond_type,
                    'Tempo': tempo,
                    'Condition_Label': cond_label,
                    'Ligament': lig,
                    'Anatomical_Group': ANATOMICAL_GROUP.get(lig, 'Other'),
                    'Resting_Length_mm': round(l0_val, 2),
                    'PCSA_Area_mm2': round(props['area'], 2),
                    'Linear_Stiffness_N': round(k_lin, 1),
                    'Viscosity_Param_eta': props['eta'],
                    'Peak_Total_Force_N': round(pk_force, 2),
                    'Mean_Total_Force_N': round(mn_force, 2),
                    'Peak_Elastic_Force_N': round(pk_elast, 2),
                    'Peak_Viscous_Force_N': round(pk_visc, 2),
                    'Viscous_Force_At_Peak_N': round(visc_at_pk, 2),
                    'Viscous_Share_At_Peak_Pct': round(visc_share, 1),
                    'Timing_Of_Peak_Force_Pct': timing_pk,
                    'Peak_Strain_Pct': round(pk_strain, 2),
                    'Mean_Strain_Pct': round(mn_strain, 2),
                    'Peak_Strain_Rate_per_s': round(pk_rate, 3)
                })
                
    df_peaks = pd.DataFrame(rows)
    print(f"Extracted {len(df_peaks)} ligament peak observations.")
    return df_peaks

def build_table4(df_peaks):
    print("\nBUILDING TABLE 4: PARTICIPANT AND COHORT PEAK METRICS")
    
    # Primary manuscript ligaments: pRCL is primary functional lateral bundle
    manuscript_ligs = ['pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'pmCAP', 'plCAP']
    df_man = df_peaks[df_peaks['Ligament'].isin(manuscript_ligs)].copy()
    
    # Also save Supplementary Table S4 for all lateral sub-bundles
    sub_bundle_ligs = ['pRCL', 'cRCL', 'aRCL', 'RCL']
    df_sub = df_peaks[df_peaks['Ligament'].isin(sub_bundle_ligs)].copy()
    sub_out = os.path.join(TABLES_DIR, "Supplementary_Table_S4_Sub_Bundles.xlsx")
    df_sub.to_excel(sub_out, index=False)
    print(f"Saved Supplementary Table S4 to {sub_out}")
    
    # Cohort Mean and SD grouped by Ligament x Condition
    cohort_summary = df_man.groupby(['Ligament', 'Condition', 'Condition_Type', 'Tempo']).agg({
        'Peak_Total_Force_N': ['mean', 'std'],
        'Peak_Elastic_Force_N': ['mean', 'std'],
        'Peak_Viscous_Force_N': ['mean', 'std'],
        'Viscous_Share_At_Peak_Pct': ['mean', 'std'],
        'Peak_Strain_Pct': ['mean', 'std'],
        'Peak_Strain_Rate_per_s': ['mean', 'std']
    }).reset_index()
    
    # Flatten multi-level column index
    cohort_summary.columns = [
        '_'.join(c).strip('_') for c in cohort_summary.columns.values
    ]
    
    # Create formatted strings (Mean ± SD)
    cohort_formatted = pd.DataFrame({
        'Ligament': cohort_summary['Ligament'],
        'Condition_Type': cohort_summary['Condition_Type'],
        'Tempo': cohort_summary['Tempo'],
        'Condition': cohort_summary['Condition'],
        'Peak_Total_Force_N': cohort_summary.apply(lambda r: f"{r['Peak_Total_Force_N_mean']:.2f} ± {r['Peak_Total_Force_N_std']:.2f}", axis=1),
        'Peak_Elastic_Force_N': cohort_summary.apply(lambda r: f"{r['Peak_Elastic_Force_N_mean']:.2f} ± {r['Peak_Elastic_Force_N_std']:.2f}", axis=1),
        'Peak_Viscous_Force_N': cohort_summary.apply(lambda r: f"{r['Peak_Viscous_Force_N_mean']:.2f} ± {r['Peak_Viscous_Force_N_std']:.2f}", axis=1),
        'Viscous_Share_Pct': cohort_summary.apply(lambda r: f"{r['Viscous_Share_At_Peak_Pct_mean']:.1f} ± {r['Viscous_Share_At_Peak_Pct_std']:.1f}%", axis=1),
        'Peak_Strain_Pct': cohort_summary.apply(lambda r: f"{r['Peak_Strain_Pct_mean']:.2f} ± {r['Peak_Strain_Pct_std']:.2f}%", axis=1),
        'Peak_Strain_Rate_s-1': cohort_summary.apply(lambda r: f"{r['Peak_Strain_Rate_per_s_mean']:.3f} ± {r['Peak_Strain_Rate_per_s_std']:.3f}", axis=1)
    })
    
    # Export to Excel
    with pd.ExcelWriter(OUT_FILE) as writer:
        cohort_formatted.to_excel(writer, sheet_name="Table4_Cohort_Mean_SD", index=False)
        cohort_summary.to_excel(writer, sheet_name="Cohort_Raw_Numeric", index=False)
        df_peaks.to_excel(writer, sheet_name="All_Participant_Peaks", index=False)
        
    print(f"Table 4 successfully saved to: {OUT_FILE}")
    return cohort_formatted, df_peaks

def main():
    df_peaks = extract_participant_peaks()
    if not df_peaks.empty:
        build_table4(df_peaks)

if __name__ == "__main__":
    main()
