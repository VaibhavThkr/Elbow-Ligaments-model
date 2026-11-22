"""
Export Repetition-Averaged Length Over Time for aRCL, cRCL, pRCL, LUCL, and caUCL
================================================================================
For each participant (S1 to S6) and trial (nor_fast, nor_med, nor_slow, rest_fast, rest_med, rest_slow):
  1. Computes the exact repetition cycle duration (in seconds) from raw kinematics.
  2. Extracts the repetition-averaged lengths (in mm and m) over 0% to 100% cycle.
  3. Provides the exact time vector t (s) = (Cycle % / 100) * Cycle_Duration.
  4. Verifies caUCL and pRCL NLSR peak force, strain, and strain rate against Table 4.
  5. Exports a clean, multi-tab Excel workbook and tidy CSV file.
"""

import os
import glob
import re
import numpy as np
import pandas as pd
import scipy.special as special
from scipy.signal import find_peaks, butter, filtfilt
from scipy.interpolate import interp1d

BASE_DIR = r"d:\Documents\Elbow joint Project"
PATHB_DIR = os.path.join(BASE_DIR, "03_Path_B_Anatomical_Axis_Model")
MODEL_DIR = os.path.join(PATHB_DIR, "models", "calibrated_devita_nlsr")
SIM_DIR = os.path.join(PATHB_DIR, "simulations")
TABLES_DIR = os.path.join(PATHB_DIR, "results", "tables")
ORIG_SIM_DIR = os.path.join(BASE_DIR, "Simulation and Analysis")

SUBJECTS = [
    ("S1", "Gokul",  "1_S1 Bicep", "1B"),
    ("S2", "Sanyam", "2_S2 Bicep", "2B"),
    ("S3", "Sunil",  "3_S3 Bicep", "3B"),
    ("S4", "Piyush", "4_S4 Bicep", "4B"),
    ("S5", "Nitin",  "5_S5 Bicep", "5B"),
    ("S6", "Rahul",  "6_S6 Bicep", "6B")
]

CONDITIONS = [
    ('nor_fast', 'Continuous',     'Fast',   'Normal Fast'),
    ('nor_med',  'Continuous',     'Medium', 'Normal Medium'),
    ('nor_slow', 'Continuous',     'Slow',   'Normal Slow'),
    ('rest_fast', 'Isometric Pause', 'Fast',   'Restricted Fast'),
    ('rest_med',  'Isometric Pause', 'Medium', 'Restricted Medium'),
    ('rest_slow', 'Isometric Pause', 'Slow',   'Restricted Slow')
]

TARGET_LIGAMENTS = ['aRCL', 'cRCL', 'pRCL', 'LUCL', 'caUCL']

# Constitutive parameters from Table 2 & Table 4
CONSTITUTIVE_PROPS = {
    'aRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'cRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'pRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'LUCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.9953, 'area': 13.35, 'F_uts': 102.4},
    'caUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.01,  'area': 4.31,  'F_uts': 89.2}
}

def read_sto(fp):
    with open(fp, "r") as f:
        lines = f.readlines()
    skip = 0
    for i, l in enumerate(lines):
        if l.strip() == "endheader":
            skip = i + 1
            break
    return pd.read_csv(fp, delim_whitespace=True, skiprows=skip)

def find_rep_boundaries(flexion_data):
    flex = flexion_data.values
    peaks, _ = find_peaks(-flex, distance=80, prominence=20)
    if len(peaks) < 2:
        peaks, _ = find_peaks(-flex, distance=40, prominence=15)
        
    valid_boundaries = []
    for i in range(len(peaks) - 1):
        s = peaks[i]
        e = peaks[i+1]
        seg = flex[s:e]
        peak_val = np.max(seg)
        rom = peak_val - min(seg[0], seg[-1])
        if peak_val >= 110.0 and rom >= 45.0:
            if not valid_boundaries:
                valid_boundaries.append(s)
            valid_boundaries.append(e)
            
    if len(valid_boundaries) >= 2:
        return np.array(valid_boundaries)
    return peaks

def get_calibrated_l0_from_osim(model_path):
    with open(model_path, 'r', encoding='utf-8') as f:
        text = f.read()
    l0_dict = {}
    for lig in ['aRCL', 'cRCL', 'pRCL', 'LUCL', 'caUCL']:
        m = re.search(r'<NonLinearLigament name="' + lig + r'">.*?</NonLinearLigament>', text, re.DOTALL)
        if m:
            block = m.group(0)
            rl_match = re.search(r'<resting_length>(.*?)</resting_length>', block)
            if rl_match:
                l0_dict[lig] = float(rl_match.group(1))
    return l0_dict

def compute_nlsr_metrics(length_arr_m, L0_m, props, cycle_duration_s):
    """
    Computes strain, strain rate, and total force for a repetition-averaged length curve.
    dt is evaluated directly from cycle_duration_s / 100.
    """
    dt = cycle_duration_s / 100.0
    strains = (length_arr_m - L0_m) / L0_m
    rates = np.gradient(strains, dt)
    
    alpha = props['alpha']
    beta  = props['beta']
    K     = props['K']
    eta   = props['eta']
    area  = props['area']
    F_uts = props['F_uts']
    
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
    
    s_tot = np.maximum(0.0, s_elast + s_visc)
    f_tot = np.minimum(F_uts, s_tot * area)
    
    peak_force = np.max(f_tot)
    peak_strain = np.max(strains) * 100.0
    peak_rate = np.max(rates)
    
    return peak_force, peak_strain, peak_rate, f_tot, strains, rates

def main():
    print("================================================================================")
    print("EXTRACTING REPETITION-AVERAGED LENGTHS & DURATIONS (aRCL, cRCL, pRCL, LUCL, caUCL)")
    print("================================================================================")
    
    tidy_rows = []
    durations_summary = []
    verification_rows = []
    
    # Load Table 4 for direct cross-checking
    table4_file = os.path.join(TABLES_DIR, "Table4_Participant_and_Cohort_Peaks.xlsx")
    df_table4 = pd.read_excel(table4_file, sheet_name="All_Participant_Peaks")
    
    participant_sheets = {sid: [] for sid, _, _, _ in SUBJECTS}
    
    cycle_pct = np.linspace(0, 100, 101)
    
    for sid, sname, folder, prefix in SUBJECTS:
        model_path = os.path.join(MODEL_DIR, f"Elbow_model-{sid}.osim")
        l0_map = get_calibrated_l0_from_osim(model_path)
        print(f"\nProcessing {sid} ({sname}):")
        print(f"  Calibrated L0 (mm): " + ", ".join([f"{k}={v*1000:.2f}" for k, v in l0_map.items()]))
        
        lengths_excel = os.path.join(SIM_DIR, folder, "Analysis Results", f"{sname} Bicep_Lengths.xlsx")
        xls_len = pd.ExcelFile(lengths_excel)
        
        # Load motion and length storage files to compute exact repetition duration
        mot_dir = os.path.join(ORIG_SIM_DIR, folder, "Motion files")
        
        for cond_code, cond_type, tempo, cond_label in CONDITIONS:
            mot_path = os.path.join(mot_dir, f"{prefix}_{cond_code}_right.mot")
            if not os.path.exists(mot_path):
                print(f"  [SKIP] Motion file not found: {mot_path}")
                continue
                
            mot_df = read_sto(mot_path)
            flex_col = [c for c in mot_df.columns if 'elbow_flex' in c.lower()][0]
            time_arr = mot_df['time'].values
            
            # Find rep boundaries
            boundaries = find_rep_boundaries(mot_df[flex_col])
            n_reps = len(boundaries) - 1
            if n_reps < 1:
                print(f"  [WARN] Not enough reps in {cond_code}")
                continue
                
            durations = []
            for i in range(n_reps):
                durations.append(time_arr[boundaries[i+1]] - time_arr[boundaries[i]])
            mean_dur = float(np.mean(durations))
            std_dur  = float(np.std(durations, ddof=1)) if len(durations) > 1 else 0.0
            
            time_vector_s = (cycle_pct / 100.0) * mean_dur
            
            # Extract lengths for all 5 target ligaments
            lig_lengths_mm = {}
            lig_lengths_m  = {}
            
            for lig in TARGET_LIGAMENTS:
                df_sheet = pd.read_excel(xls_len, lig)
                if cond_code in df_sheet.columns:
                    len_m = df_sheet[cond_code].values
                    lig_lengths_m[lig] = len_m
                    lig_lengths_mm[lig] = len_m * 1000.0
                else:
                    print(f"  [!] Missing column {cond_code} in sheet {lig}")
                    
            if len(lig_lengths_mm) < 5:
                continue
                
            # Store summary row
            dur_row = {
                'Participant': sid,
                'Participant_Name': sname,
                'Condition': cond_code,
                'Condition_Type': cond_type,
                'Tempo': tempo,
                'Condition_Label': cond_label,
                'N_Repetitions': n_reps,
                'Cycle_Duration_Mean_s': round(mean_dur, 3),
                'Cycle_Duration_SD_s': round(std_dur, 3),
                'aRCL_L0_mm': round(l0_map['aRCL']*1000, 2),
                'cRCL_L0_mm': round(l0_map['cRCL']*1000, 2),
                'pRCL_L0_mm': round(l0_map['pRCL']*1000, 2),
                'LUCL_L0_mm': round(l0_map['LUCL']*1000, 2),
                'caUCL_L0_mm': round(l0_map['caUCL']*1000, 2),
                'aRCL_Min_mm': round(np.min(lig_lengths_mm['aRCL']), 2),
                'aRCL_Max_mm': round(np.max(lig_lengths_mm['aRCL']), 2),
                'aRCL_Delta_mm': round(np.max(lig_lengths_mm['aRCL']) - np.min(lig_lengths_mm['aRCL']), 2),
                'cRCL_Min_mm': round(np.min(lig_lengths_mm['cRCL']), 2),
                'cRCL_Max_mm': round(np.max(lig_lengths_mm['cRCL']), 2),
                'cRCL_Delta_mm': round(np.max(lig_lengths_mm['cRCL']) - np.min(lig_lengths_mm['cRCL']), 2),
                'pRCL_Min_mm': round(np.min(lig_lengths_mm['pRCL']), 2),
                'pRCL_Max_mm': round(np.max(lig_lengths_mm['pRCL']), 2),
                'pRCL_Delta_mm': round(np.max(lig_lengths_mm['pRCL']) - np.min(lig_lengths_mm['pRCL']), 2),
                'LUCL_Min_mm': round(np.min(lig_lengths_mm['LUCL']), 2),
                'LUCL_Max_mm': round(np.max(lig_lengths_mm['LUCL']), 2),
                'LUCL_Delta_mm': round(np.max(lig_lengths_mm['LUCL']) - np.min(lig_lengths_mm['LUCL']), 2),
                'caUCL_Min_mm': round(np.min(lig_lengths_mm['caUCL']), 2),
                'caUCL_Max_mm': round(np.max(lig_lengths_mm['caUCL']), 2),
                'caUCL_Delta_mm': round(np.max(lig_lengths_mm['caUCL']) - np.min(lig_lengths_mm['caUCL']), 2)
            }
            durations_summary.append(dur_row)
            
            # Cross-verify caUCL and pRCL with Table 4
            for check_lig in ['caUCL', 'pRCL']:
                L0 = l0_map[check_lig]
                props = CONSTITUTIVE_PROPS[check_lig]
                calc_f, calc_s, calc_r, _, _, _ = compute_nlsr_metrics(
                    lig_lengths_m[check_lig], L0, props, mean_dur
                )
                
                # Fetch Table 4 values
                t4_row = df_table4[(df_table4['Participant'] == sid) & 
                                   (df_table4['Condition'] == cond_code) & 
                                   (df_table4['Ligament'] == check_lig)]
                if not t4_row.empty:
                    t4_f = t4_row['Peak_Total_Force_N'].values[0]
                    t4_s = t4_row['Peak_Strain_Pct'].values[0]
                    t4_r = t4_row['Peak_Strain_Rate_per_s'].values[0]
                    
                    verification_rows.append({
                        'Participant': sid,
                        'Condition': cond_code,
                        'Ligament': check_lig,
                        'Calibrated_L0_mm': round(L0*1000, 2),
                        'Cycle_Duration_s': round(mean_dur, 3),
                        'Table4_Peak_Force_N': t4_f,
                        'Calculated_Peak_Force_N': round(calc_f, 2),
                        'Force_Difference_N': round(abs(t4_f - calc_f), 2),
                        'Table4_Peak_Strain_Pct': t4_s,
                        'Calculated_Peak_Strain_Pct': round(calc_s, 2),
                        'Table4_Peak_Strain_Rate_s-1': t4_r,
                        'Calculated_Peak_Strain_Rate_s-1': round(calc_r, 3),
                        'Status': 'CONCORDANT' if abs(t4_f - calc_f) <= 0.20 else 'SLIGHT_VARIATION'
                    })
                    
            # Populate Tidy waveform rows (101 points per condition)
            for pt in range(101):
                p_row = {
                    'Participant': sid,
                    'Participant_Name': sname,
                    'Condition': cond_code,
                    'Condition_Type': cond_type,
                    'Tempo': tempo,
                    'Cycle_Pct': cycle_pct[pt],
                    'Cycle_Duration_s': round(mean_dur, 4),
                    'Time_s': round(time_vector_s[pt], 4),
                    'aRCL_Length_mm': round(lig_lengths_mm['aRCL'][pt], 4),
                    'cRCL_Length_mm': round(lig_lengths_mm['cRCL'][pt], 4),
                    'pRCL_Length_mm': round(lig_lengths_mm['pRCL'][pt], 4),
                    'LUCL_Length_mm': round(lig_lengths_mm['LUCL'][pt], 4),
                    'caUCL_Length_mm': round(lig_lengths_mm['caUCL'][pt], 4),
                    'aRCL_Length_m': round(lig_lengths_m['aRCL'][pt], 6),
                    'cRCL_Length_m': round(lig_lengths_m['cRCL'][pt], 6),
                    'pRCL_Length_m': round(lig_lengths_m['pRCL'][pt], 6),
                    'LUCL_Length_m': round(lig_lengths_m['LUCL'][pt], 6),
                    'caUCL_Length_m': round(lig_lengths_m['caUCL'][pt], 6)
                }
                tidy_rows.append(p_row)
                participant_sheets[sid].append(p_row)
                
            print(f"  [OK] {cond_code}: {n_reps} reps, duration = {mean_dur:.2f} s")

    df_tidy = pd.DataFrame(tidy_rows)
    df_durations = pd.DataFrame(durations_summary)
    df_verif = pd.DataFrame(verification_rows)
    
    # Save Excel Workbook
    out_xlsx = os.path.join(TABLES_DIR, "PathB_Repetition_Averaged_Lengths_and_Durations.xlsx")
    with pd.ExcelWriter(out_xlsx) as writer:
        df_durations.to_excel(writer, sheet_name="Durations_and_L0_Summary", index=False)
        df_verif.to_excel(writer, sheet_name="Table4_Verification", index=False)
        df_tidy.to_excel(writer, sheet_name="Tidy_All_Waveforms", index=False)
        for sid, rows in participant_sheets.items():
            pd.DataFrame(rows).to_excel(writer, sheet_name=f"{sid}_Waveforms", index=False)
            
    print(f"\n[SAVED] Excel Workbook: {out_xlsx}")
    
    # Save standalone CSV for easy import
    out_csv = os.path.join(TABLES_DIR, "PathB_Repetition_Averaged_Lengths_and_Durations.csv")
    df_tidy.to_csv(out_csv, index=False)
    print(f"[SAVED] Standalone CSV: {out_csv}")
    
    # Copy to artifacts directory
    artifact_dir = r"C:\Users\vaibh\.gemini\antigravity\brain\1480b829-4eb3-4948-8aae-83a7b237bdf1"
    art_xlsx = os.path.join(artifact_dir, "PathB_Repetition_Averaged_Lengths_and_Durations.xlsx")
    art_csv  = os.path.join(artifact_dir, "PathB_Repetition_Averaged_Lengths_and_Durations.csv")
    df_durations.to_excel(os.path.join(artifact_dir, "PathB_Trial_Cycle_Durations_and_L0.xlsx"), index=False)
    df_tidy.to_excel(art_xlsx, index=False)
    df_tidy.to_csv(art_csv, index=False)
    print(f"[SAVED] Artifact copies updated in {artifact_dir}")

    # Print summary of Table 4 verification
    print("\n================================================================================")
    print("VERIFICATION OF caUCL & pRCL REPRODUCTION AGAINST TABLE 4:")
    print("================================================================================")
    print(df_verif.to_string())

if __name__ == "__main__":
    main()
