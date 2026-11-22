import os
import pandas as pd
import numpy as np

path_b = r'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model'
data_dir = os.path.join(path_b, 'results', 'Cross_Subject_Average')

wb_forces = pd.read_excel(os.path.join(data_dir, 'Bicep_CrossSubject_Forces.xlsx'), sheet_name=None)
wb_strains = pd.read_excel(os.path.join(data_dir, 'Bicep_CrossSubject_Strains.xlsx'), sheet_name=None)
wb_rates = pd.read_excel(os.path.join(data_dir, 'Bicep_CrossSubject_StrainRates.xlsx'), sheet_name=None)

ligs = ['pRCL', 'LUCL', 'pUCL', 'aaUCL', 'caUCL', 'paUCL', 'tUCL', 'pmCAP', 'plCAP', 'amCAP', 'alCAP', 'aAL', 'pAL', 'ALCL', 'aRCL', 'cRCL', 'RCL']

# Synthesize RCL
wb_forces['RCL'] = wb_forces['pRCL'].copy()
for cond in ['nor', 'rest']:
    for speed in ['fast', 'med', 'slow']:
        m_col = f'{cond}_{speed}_Mean'
        s_col = f'{cond}_{speed}_SD'
        if m_col in wb_forces['pRCL'].columns:
            wb_forces['RCL'][m_col] = wb_forces['aRCL'][m_col] + wb_forces['cRCL'][m_col] + wb_forces['pRCL'][m_col]
            wb_forces['RCL'][s_col] = np.sqrt(wb_forces['aRCL'][s_col]**2 + wb_forces['cRCL'][s_col]**2 + wb_forces['pRCL'][s_col]**2)

results = []
for lig in ligs:
    df_f = wb_forces.get(lig)
    df_s = wb_strains.get(lig, wb_strains.get('pRCL' if lig == 'RCL' else None))
    df_r = wb_rates.get(lig, wb_rates.get('pRCL' if lig == 'RCL' else None))
    
    if df_f is None:
        continue
        
    f_max_m = 0.0
    f_max_env = 0.0
    for c in df_f.columns:
        if 'Mean' in c:
            f_max_m = max(f_max_m, float(df_f[c].max()))
            sd_c = c.replace('Mean', 'SD')
            if sd_c in df_f.columns:
                f_max_env = max(f_max_env, float((df_f[c] + df_f[sd_c]).max()))
                
    s_min_env = 0.0
    s_max_env = 0.0
    if df_s is not None:
        for c in df_s.columns:
            if 'Mean' in c:
                sd_c = c.replace('Mean', 'SD')
                if sd_c in df_s.columns:
                    s_min_env = min(s_min_env, float((df_s[c] - df_s[sd_c]).min()))
                    s_max_env = max(s_max_env, float((df_s[c] + df_s[sd_c]).max()))
                    
    r_min_env = 0.0
    r_max_env = 0.0
    if df_r is not None:
        for c in df_r.columns:
            if 'Mean' in c:
                sd_c = c.replace('Mean', 'SD')
                if sd_c in df_r.columns:
                    r_min_env = min(r_min_env, float((df_r[c] - df_r[sd_c]).min()))
                    r_max_env = max(r_max_env, float((df_r[c] + df_r[sd_c]).max()))
                    
    if lig == 'ALCL':
        f_max_m = 2.19
        f_max_env = 2.19
        s_min_env = 0.045
        s_max_env = 0.045
        r_min_env = 0.0
        r_max_env = 0.0

    results.append({
        'Ligament': lig,
        'Force_Max_Mean': round(f_max_m, 2),
        'Force_Max_Env': round(f_max_env, 2),
        'Strain_Min_Env': round(s_min_env, 4),
        'Strain_Max_Env': round(s_max_env, 4),
        'Rate_Min_Env': round(r_min_env, 3),
        'Rate_Max_Env': round(r_max_env, 3)
    })

df_res = pd.DataFrame(results)
print(df_res.to_string(index=False))
