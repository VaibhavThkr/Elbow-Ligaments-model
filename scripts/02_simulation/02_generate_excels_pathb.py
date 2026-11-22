"""
Phase 4: Excel Data Generation for Path B Anatomical Axis Model
==============================================================
Extracts repetition-averaged, time-normalized (0% to 100% cycle) data from the
generated .sto simulation files across all 35 trials and 6 subjects.

Integrates Option A: 3-Fiber RCL Fan (aRCL, cRCL, pRCL) + LUCL + Medial + Capsular structures.

Generates 4 Excel workbooks per subject:
  1. <Subject>_Forces.xlsx
  2. <Subject>_Lengths.xlsx
  3. <Subject>_Strains.xlsx
  4. <Subject>_StrainRates.xlsx

And 4 Cross-Subject Grand Average workbooks:
  1. Bicep_CrossSubject_Forces.xlsx
  2. Bicep_CrossSubject_Lengths.xlsx
  3. Bicep_CrossSubject_Strains.xlsx
  4. Bicep_CrossSubject_StrainRates.xlsx
"""

import os
import glob
import re
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, butter, filtfilt
from scipy.interpolate import interp1d
from scipy import special

BASE_DIR = r"d:\Documents\Elbow joint Project"
PATHB_DIR = os.path.join(BASE_DIR, "03_Path_B_Anatomical_Axis_Model")
MODEL_DIR = os.path.join(PATHB_DIR, "models", "calibrated_devita_nlsr")
SIM_DIR = os.path.join(PATHB_DIR, "simulations")
RESULTS_DIR = os.path.join(PATHB_DIR, "results")
CROSS_DIR = os.path.join(RESULTS_DIR, "Cross_Subject_Average")
os.makedirs(CROSS_DIR, exist_ok=True)

# Full ligament set including 3-fiber RCL fan
LIGAMENTS = [
    'RCL', 'aRCL', 'cRCL', 'pRCL', 'LUCL',
    'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'tUCL',
    'pmCAP', 'plCAP', 'amCAP', 'alCAP',
    'aAL', 'pAL', 'ALCL', 'QL', 'OC'
]

SUBJECTS = [
    ("S1", "Gokul",  "1_S1 Bicep", "1B"),
    ("S2", "Sanyam", "2_S2 Bicep", "2B"),
    ("S3", "Sunil",  "3_S3 Bicep", "3B"),
    ("S4", "Piyush", "4_S4 Bicep", "4B"),
    ("S5", "Nitin",  "5_S5 Bicep", "5B"),
    ("S6", "Rahul",  "6_S6 Bicep", "6B")
]

FILTER_CUTOFF = 15.0  # Hz
FILTER_ORDER = 4

# Path B Calibrated Constitutive Properties
LIGAMENT_PROPS = {
    # 3-Fiber RCL Fan
    'aRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'cRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'pRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67},
    'RCL':   {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 17.01, 'F_uts': 110.0},
    # Calibrated LUCL
    'LUCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.9953, 'area': 13.35, 'F_uts': 102.4},
    # Medial UCL
    'aaUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.81,  'area': 4.31,  'F_uts': 89.2},
    'caUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.01,  'area': 4.31,  'F_uts': 89.2},
    'paUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 7.35,  'area': 4.31,  'F_uts': 89.2},
    'pUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 8.90,  'F_uts': 160.0},
    'tUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 2.23,  'F_uts': 160.0},
    # Posterior & Anterior Capsule
    'pmCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.62,  'area': 4.39,  'F_uts': 126.3},
    'plCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.96,  'area': 4.39,  'F_uts': 126.3},
    'amCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.62,  'area': 4.39,  'F_uts': 126.3},
    'alCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.96,  'area': 4.39,  'F_uts': 126.3},
    # Annular & Accessory
    'aAL':   {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.94,  'area': 6.48,  'F_uts': 102.95},
    'pAL':   {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 9.87,  'area': 6.48,  'F_uts': 102.95},
    'ALCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 2.44,  'area': 4.25,  'F_uts': 65.0},
    'QL':    {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 0.00,  'area': 0.00,  'F_uts': 0.0},
    'OC':    {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 0.00,  'area': 0.00,  'F_uts': 0.0}
}

def compute_nlsr_force(eps, epsdot, props):
    """Computes tension-only rectified physiological NLSR force (F >= 0)."""
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
    return F_tot

def read_opensim_storage(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    skip_lines = 0
    for i, line in enumerate(lines):
        if line.strip() == 'endheader':
            skip_lines = i + 1
            break
    df = pd.read_csv(filepath, sep='\t', skiprows=skip_lines)
    df.columns = [col.strip() for col in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def get_sampling_frequency(time_array):
    dt = np.mean(np.diff(time_array))
    if dt == 0 or np.isnan(dt):
        return 100.0
    return 1.0 / dt

def get_ligament_column(df, lig_name):
    low = lig_name.lower()
    for col in df.columns:
        c_low = col.lower()
        if c_low == low or c_low == f"{low}_m" or c_low.startswith(f"{low}_"):
            return col
    return None

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
    for lig in LIGAMENTS:
        m = re.search(r'<NonLinearLigament name="' + lig + r'">.*?</NonLinearLigament>', text, re.DOTALL)
        if m:
            block = m.group(0)
            rl_match = re.search(r'<resting_length>(.*?)</resting_length>', block)
            if rl_match:
                l0_dict[lig] = float(rl_match.group(1))
    # Fallbacks for Path B sub-fibers if not directly named in NonLinearLigament
    if 'aRCL' not in l0_dict: l0_dict['aRCL'] = 0.02000
    if 'cRCL' not in l0_dict: l0_dict['cRCL'] = 0.01973
    if 'pRCL' not in l0_dict: l0_dict['pRCL'] = 0.02030
    if 'LUCL' not in l0_dict: l0_dict['LUCL'] = 0.05250
    if 'RCL'  not in l0_dict: l0_dict['RCL']  = 0.01981
    return l0_dict

def process_subject(sid, name, bicep_folder, prefix):
    print(f"\n=======================================================")
    print(f"PROCESSING PATH B DATA FOR: {sid} ({name.upper()})")
    print(f"=======================================================")
    
    analysis_dir = os.path.join(SIM_DIR, bicep_folder, "Analysis Results")
    motion_dir = os.path.join(BASE_DIR, "Simulation and Analysis", bicep_folder, "Motion files")
    model_path = os.path.join(MODEL_DIR, f"Elbow_model-{sid}.osim")
    
    if not os.path.exists(analysis_dir) or not os.path.exists(motion_dir):
        print(f"  [ERROR] Missing analysis or motion directory for {sid}")
        return
        
    l0_map = get_calibrated_l0_from_osim(model_path)
    print(f"  Extracted {len(l0_map)} calibrated resting lengths from {os.path.basename(model_path)}")
    
    mot_files = sorted(glob.glob(os.path.join(motion_dir, "*.mot")))
    conditions = [os.path.splitext(os.path.basename(f))[0] for f in mot_files]
    print(f"  Found {len(conditions)} conditions: {[c.replace(prefix+'_', '') for c in conditions]}")
    
    grand_means = {
        'Force':      {lig: {} for lig in LIGAMENTS},
        'Length':     {lig: {} for lig in LIGAMENTS},
        'Strain':     {lig: {} for lig in LIGAMENTS},
        'StrainRate': {lig: {} for lig in LIGAMENTS},
    }
    
    for condition in conditions:
        mot_path = os.path.join(motion_dir, f"{condition}.mot")
        trial_dir = os.path.join(analysis_dir, condition)
        
        force_files = glob.glob(os.path.join(trial_dir, "*_ForceReporter_forces.sto"))
        length_files = glob.glob(os.path.join(trial_dir, "*_MuscleAnalysis_Length.sto"))
        
        if not os.path.exists(mot_path) or not force_files or not length_files:
            print(f"    [SKIP] Missing .sto or .mot for {condition}")
            continue
            
        mot_df = read_opensim_storage(mot_path)
        force_df = read_opensim_storage(force_files[0])
        length_df = read_opensim_storage(length_files[0])
        
        # Identify flexion column
        flex_col = None
        for col in mot_df.columns:
            c_low = col.lower()
            if 'elbow_flex_r' in c_low and 'right' in condition.lower():
                flex_col = col; break
            elif 'elbow_flex_l' in c_low and 'left' in condition.lower():
                flex_col = col; break
            elif 'flexion' in c_low:
                flex_col = col; break
        if flex_col is None:
            for col in mot_df.columns:
                if 'elbow_flex' in col.lower():
                    flex_col = col; break
                    
        if flex_col is None:
            print(f"    [!] Flexion column not found in {mot_path}")
            continue
            
        minima_indices = find_rep_boundaries(mot_df[flex_col])
        if len(minima_indices) < 2:
            print(f"    [!] Not enough repetitions detected for {condition} ({len(minima_indices)} found)")
            continue
            
        time_array = length_df['time'].values
        fs = get_sampling_frequency(time_array)
        b, a = butter(FILTER_ORDER, FILTER_CUTOFF, fs=fs, btype='lowpass')
        dt = 1.0 / fs

        # Pre-process all individual ligament waveforms
        computed_signals = {}
        for lig in LIGAMENTS:
            if lig == 'RCL': continue # Processed as sum of fan fibers
            len_col = get_ligament_column(length_df, lig)
            if not len_col:
                continue
            raw_len = length_df[len_col].values
            raw_len = np.nan_to_num(raw_len, nan=0.0)
            filt_len = filtfilt(b, a, raw_len)
            L0 = l0_map.get(lig, 0.02)
            filt_strain = (filt_len - L0) / L0
            if lig == 'ALCL':
                filt_strain = np.zeros_like(filt_len)
                s_rate = np.zeros_like(filt_len)
            else:
                s_rate = np.gradient(filt_strain, dt)
            props = LIGAMENT_PROPS.get(lig, LIGAMENT_PROPS['RCL'])
            filt_force = compute_nlsr_force(filt_strain, s_rate, props)
            computed_signals[lig] = {
                'length': filt_len,
                'strain': filt_strain,
                'strain_rate': s_rate,
                'force': filt_force
            }

        # Synthesize combined RCL from 3-fiber fan
        if 'aRCL' in computed_signals and 'cRCL' in computed_signals and 'pRCL' in computed_signals:
            rcl_tot_force = (computed_signals['aRCL']['force'] + 
                             computed_signals['cRCL']['force'] + 
                             computed_signals['pRCL']['force'])
            rcl_mean_len = (computed_signals['aRCL']['length'] + 
                            computed_signals['cRCL']['length'] + 
                            computed_signals['pRCL']['length']) / 3.0
            L0_rcl = l0_map.get('RCL', 0.01981)
            rcl_mean_strain = (rcl_mean_len - L0_rcl) / L0_rcl
            rcl_mean_srate = np.gradient(rcl_mean_strain, dt)
            computed_signals['RCL'] = {
                'length': rcl_mean_len,
                'strain': rcl_mean_strain,
                'strain_rate': rcl_mean_srate,
                'force': rcl_tot_force
            }
        elif get_ligament_column(length_df, 'RCL'):
            # Legacy fallback if individual fan columns are not present
            len_col = get_ligament_column(length_df, 'RCL')
            raw_len = length_df[len_col].values
            filt_len = filtfilt(b, a, np.nan_to_num(raw_len))
            L0 = l0_map.get('RCL', 0.01981)
            filt_strain = (filt_len - L0) / L0
            s_rate = np.gradient(filt_strain, dt)
            filt_force = compute_nlsr_force(filt_strain, s_rate, LIGAMENT_PROPS['RCL'])
            computed_signals['RCL'] = {
                'length': filt_len, 'strain': filt_strain,
                'strain_rate': s_rate, 'force': filt_force
            }

        # Time normalization across repetitions
        for lig, sigs in computed_signals.items():
            reps_force, reps_length, reps_strain, reps_srate = [], [], [], []
            for i in range(len(minima_indices) - 1):
                s_idx = minima_indices[i]
                e_idx = minima_indices[i+1]
                max_len = len(sigs['length'])
                if e_idx > max_len: e_idx = max_len
                if s_idx >= max_len: break
                
                rep_t = time_array[s_idx:e_idx]
                if len(rep_t) < 2: continue
                
                if 'rest' in condition.lower():
                    rep_flex = mot_df[flex_col].values[s_idx:e_idx]
                    threshold = np.max(rep_flex) - 5.0
                    above_thresh = np.where(rep_flex > threshold)[0]
                    if len(above_thresh) >= 2 and above_thresh[0] > 0 and above_thresh[-1] < len(rep_flex) - 1:
                        s_hold, e_hold = above_thresh[0], above_thresh[-1]
                        t_conc = np.linspace(0, 0.33, s_hold + 1)
                        t_hold = np.linspace(0.33, 0.67, e_hold - s_hold + 1)
                        t_ecc  = np.linspace(0.67, 1.0, len(rep_flex) - e_hold)
                        rep_t_norm = np.concatenate([t_conc[:-1], t_hold[:-1], t_ecc])
                    else:
                        rep_t_norm = (rep_t - rep_t[0]) / (rep_t[-1] - rep_t[0])
                else:
                    rep_t_norm = (rep_t - rep_t[0]) / (rep_t[-1] - rep_t[0])
                    
                t_target = np.linspace(0, 1, 101)
                f_f = interp1d(rep_t_norm, sigs['force'][s_idx:e_idx], kind='linear', fill_value="extrapolate")
                f_l = interp1d(rep_t_norm, sigs['length'][s_idx:e_idx], kind='linear', fill_value="extrapolate")
                f_s = interp1d(rep_t_norm, sigs['strain'][s_idx:e_idx], kind='linear', fill_value="extrapolate")
                f_r = interp1d(rep_t_norm, sigs['strain_rate'][s_idx:e_idx], kind='linear', fill_value="extrapolate")
                
                reps_force.append(f_f(t_target))
                reps_length.append(f_l(t_target))
                reps_strain.append(f_s(t_target))
                reps_srate.append(f_r(t_target))
                
            if reps_force:
                grand_means['Force'][lig][condition]      = np.mean(reps_force, axis=0)
                grand_means['Length'][lig][condition]     = np.mean(reps_length, axis=0)
                grand_means['Strain'][lig][condition]     = np.mean(reps_strain, axis=0)
                grand_means['StrainRate'][lig][condition] = np.mean(reps_srate, axis=0)
                
        print(f"    [OK] Processed {condition} ({len(minima_indices)-1} reps)")

    export_subject_excels(analysis_dir, sid, name, grand_means, conditions, prefix)

def export_subject_excels(analysis_dir, sid, name, grand_means, conditions, prefix):
    file_map = {
        'Force':      f"{name} Bicep_Forces.xlsx",
        'Length':     f"{name} Bicep_Lengths.xlsx",
        'Strain':     f"{name} Bicep_Strains.xlsx",
        'StrainRate': f"{name} Bicep_StrainRates.xlsx"
    }
    
    cycle_pct = np.linspace(0, 100, 101)
    for var, filename in file_map.items():
        out_path = os.path.join(analysis_dir, filename)
        has_data = any(grand_means[var][l] for l in LIGAMENTS if l in grand_means[var])
        if not has_data:
            print(f"  [!] No data found for {var}, skipping {filename}")
            continue
            
        with pd.ExcelWriter(out_path) as writer:
            for lig in LIGAMENTS:
                if lig not in grand_means[var] or not grand_means[var][lig]:
                    continue
                df_data = {'Cycle %': cycle_pct}
                for condition in conditions:
                    if condition in grand_means[var][lig]:
                        col_name = condition.replace(prefix + "_", "").replace("_right", "")
                        df_data[col_name] = grand_means[var][lig][condition]
                df = pd.DataFrame(df_data)
                df.to_excel(writer, sheet_name=lig, index=False)
                
        print(f"  [SAVED] {filename} -> {analysis_dir}")

def export_cross_subject_excels():
    print("\n==================================================================")
    print("COMPUTING PATH B CROSS-SUBJECT GRAND AVERAGES")
    print("==================================================================")
    
    variables = ['Forces', 'Lengths', 'Strains', 'StrainRates']
    var_keys  = ['Force', 'Length', 'Strain', 'StrainRate']
    speeds = ['fast', 'med', 'slow']
    conditions = ['nor', 'rest']
    cycle_pct = np.linspace(0, 100, 101)
    
    for var in variables:
        subject_dfs = {}
        for sid, name, bicep_folder, prefix in SUBJECTS:
            fpath = os.path.join(SIM_DIR, bicep_folder, "Analysis Results", f"{name} Bicep_{var}.xlsx")
            if os.path.exists(fpath):
                subject_dfs[sid] = pd.ExcelFile(fpath)
                
        if not subject_dfs:
            print(f"  [!] No subject workbooks found for {var}")
            continue
            
        out_excel = os.path.join(CROSS_DIR, f"Bicep_CrossSubject_{var}.xlsx")
        with pd.ExcelWriter(out_excel) as writer:
            for lig in LIGAMENTS:
                lig_cols = {'Cycle %': cycle_pct}
                for cond in conditions:
                    for spd in speeds:
                        col_key = f"{cond}_{spd}"
                        stacked = []
                        for sid, xls in subject_dfs.items():
                            if lig in xls.sheet_names:
                                s_df = pd.read_excel(xls, sheet_name=lig)
                                if col_key in s_df.columns:
                                    vals = s_df[col_key].values
                                    if len(vals) == 101:
                                        stacked.append(vals)
                        if stacked:
                            arr = np.array(stacked)
                            lig_cols[f"{cond}_{spd}_Mean"] = np.nanmean(arr, axis=0)
                            lig_cols[f"{cond}_{spd}_SD"]   = np.nanstd(arr, axis=0, ddof=1)
                if len(lig_cols) > 1:
                    df_cross = pd.DataFrame(lig_cols)
                    df_cross.to_excel(writer, sheet_name=lig, index=False)
        print(f"  [SAVED] Cross-Subject Workbook: {out_excel}")

def main():
    print("==================================================================")
    print("STARTING PATH B PARTICIPANT EXCEL GENERATION (ALL 6 PARTICIPANTS)")
    print("==================================================================")
    for sid, name, bicep_folder, prefix in SUBJECTS:
        process_subject(sid, name, bicep_folder, prefix)
    export_cross_subject_excels()
    print("\n==================================================================")
    print("PATH B EXCEL GENERATION PIPELINE COMPLETE!")
    print("==================================================================")

if __name__ == "__main__":
    main()
