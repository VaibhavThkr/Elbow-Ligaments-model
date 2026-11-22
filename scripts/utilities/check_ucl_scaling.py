import re, glob, os
import pandas as pd

models = sorted(glob.glob(r'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model\models\calibrated_devita_nlsr\Elbow_model-S*.osim'))
print(f"Found {len(models)} Path B models.")

for mpath in models:
    fname = os.path.basename(mpath)
    with open(mpath, 'r', encoding='utf-8') as f:
        txt = f.read()
    
    print(f"\n--- {fname} ---")
    for lig in ['aaUCL', 'caUCL', 'paUCL', 'pUCL']:
        m = re.search(r'<NonLinearLigament name="' + lig + r'">(.*?)</NonLinearLigament>', txt, re.DOTALL)
        if m:
            pts = re.findall(r'<PathPoint name="([^"]+)">.*?<location>([^<]+)</location>', m.group(1), re.DOTALL)
            print(f"  {lig}: {pts}")

# Also check original models in Simulation and Analysis
orig_models = sorted(glob.glob(r'd:\Documents\Elbow joint Project\Simulation and Analysis\*\Elbow_model-*.osim'))
print(f"\nFound {len(orig_models)} original models.")
for mpath in orig_models:
    fname = os.path.basename(mpath)
    folder = os.path.basename(os.path.dirname(mpath))
    with open(mpath, 'r', encoding='utf-8') as f:
        txt = f.read()
    print(f"\n--- Original {folder} / {fname} ---")
    for lig in ['aaUCL', 'caUCL', 'paUCL', 'pUCL']:
        m = re.search(r'<NonLinearLigament name="' + lig + r'">(.*?)</NonLinearLigament>', txt, re.DOTALL)
        if m:
            pts = re.findall(r'<PathPoint name="([^"]+)">.*?<location>([^<]+)</location>', m.group(1), re.DOTALL)
            print(f"  {lig}: {pts}")
        else:
            # Check Millard or other element
            m2 = re.search(r'<Millard2012EquilibriumMuscle name="' + lig + r'_m">(.*?)</Millard2012EquilibriumMuscle>', txt, re.DOTALL)
            if m2:
                pts = re.findall(r'<PathPoint name="([^"]+)">.*?<location>([^<]+)</location>', m2.group(1), re.DOTALL)
                print(f"  {lig}_m: {pts}")
