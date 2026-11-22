import os, glob
import pandas as pd
import numpy as np

base = r'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model\simulations'
subjects = [
    ('S1', '1_S1 Bicep'),
    ('S2', '2_S2 Bicep'),
    ('S3', '3_S3 Bicep'),
    ('S4', '4_S4 Bicep'),
    ('S5', '5_S5 Bicep'),
    ('S6', '6_S6 Bicep')
]

ligs = ['cRCL_m', 'aRCL_m', 'pRCL_m', 'LUCL_m', 'aaUCL_m', 'caUCL_m', 'paUCL_m', 'pUCL_m']

for sid, sdir in subjects:
    s_path = os.path.join(base, sdir, 'Analysis Results')
    length_files = glob.glob(os.path.join(s_path, '*', '*_MuscleAnalysis_Length.sto'))
    
    stats = {lig: {'l_start': [], 'l_max': [], 'l_min': []} for lig in ligs}
    
    for lf in length_files:
        with open(lf, 'r') as f:
            for i in range(40):
                line = f.readline()
                if 'endheader' in line:
                    end_idx = i + 1
                    break
        df = pd.read_csv(lf, delim_whitespace=True, skiprows=end_idx)
        
        for lig in ligs:
            if lig in df.columns:
                stats[lig]['l_start'].append(df[lig].iloc[0])
                stats[lig]['l_max'].append(df[lig].max())
                stats[lig]['l_min'].append(df[lig].min())
                
    print(f'=== {sid} ({len(length_files)} trials) ===', flush=True)
    for lig in ligs:
        l_st = np.mean(stats[lig]['l_start']) * 1000 if stats[lig]['l_start'] else 0
        l_mx = np.mean(stats[lig]['l_max']) * 1000 if stats[lig]['l_max'] else 0
        l_mn = np.mean(stats[lig]['l_min']) * 1000 if stats[lig]['l_min'] else 0
        print(f'  {lig:8s}: L_start={l_st:6.2f} mm, L_min={l_mn:6.2f} mm, L_max={l_mx:6.2f} mm, delta={l_mx-l_mn:5.2f} mm', flush=True)
