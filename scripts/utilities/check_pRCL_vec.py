import pandas as pd
import numpy as np

df = pd.read_excel(r'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model\results\Cross_Subject_Average\Bicep_CrossSubject_Lengths.xlsx', 'pRCL')
L = df['nor_fast_Mean'].values
print(f"L min: {L.min()*1000:.3f} mm, L max: {L.max()*1000:.3f} mm, L start: {L[0]*1000:.3f} mm")
