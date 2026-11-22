import os, glob
import pandas as pd
import numpy as np

base = r'PROJECT_DIR\simulations'
subjects = ['1_S1 Bicep', '2_S2 Bicep', '3_S3 Bicep', '4_S4 Bicep', '5_S5 Bicep', '6_S6 Bicep']

for s in subjects:
    s_path = os.path.join(base, s, 'Analysis Results')
    length_files = glob.glob(os.path.join(s_path, '*', '*_MuscleAnalysis_Length.sto'))
    print(f'=== {s}: {len(length_files)} length files ===')
    if length_files:
        # Check first trial
        lf = length_files[0]
        # read header
        with open(lf, 'r') as f:
            lines = [f.readline() for _ in range(30)]
        end_idx = 0
        for i, l in enumerate(lines):
            if 'endheader' in l: end_idx = i + 1; break
        df = pd.read_csv(lf, sep=r'\s+', skiprows=end_idx)
        print('  Columns:', [c for c in df.columns if any(k in c for k in ['RCL', 'LUCL'])])
        for lig in ['cRCL_m', 'aRCL_m', 'pRCL_m', 'LUCL_m']:
            if lig in df.columns:
                print(f'    {lig}: start={df[lig].iloc[0]*1000:.3f} mm, min={df[lig].min()*1000:.3f} mm, max={df[lig].max()*1000:.3f} mm')
