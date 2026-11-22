"""
Phase 1: Update Subject-Specific Models (S1-S6)
==============================================
1. Scales lateral attachment coordinates per participant using subject marker scaling.
2. Updates biophysical parameters:
   - LUCL: K = 49.2 MPa, area = 13.35 mm^2, eta* = 4.9953 MPa*s, F_uts = 102.4 N
   - aRCL: K = 49.2 MPa, area = 5.67 mm^2,  eta* = 3.7570 MPa*s, F_uts = 36.67 N
   - cRCL: K = 49.2 MPa, area = 5.67 mm^2,  eta* = 3.7570 MPa*s, F_uts = 36.67 N
   - pRCL: K = 49.2 MPa, area = 5.67 mm^2,  eta* = 3.7570 MPa*s, F_uts = 36.67 N
   - tUCL: area = 2.23 mm^2 (restored per Table 2)
   - ALCL: area = 4.25 mm^2 (restored per Table 2)
3. Calibrates L0 per participant:
   - LUCL: L0 = L_start / 1.000 (~52.9 mm for S1, eps_0 = 0.0%)
   - aRCL: L0 = L_start / 1.005 (~20.33 mm for S1, eps_0 = 0.5%)
   - cRCL: L0 = L_start (scaled axis length, eps = 0.0% throughout)
   - pRCL: L0 = L_max / 1.0487 (~20.32 mm for S1, eps_peak = 4.87%)
4. Updates both <NonLinearLigament> and <Millard2012EquilibriumMuscle> in all 6 .osim files.
"""

import os
import re
import numpy as np
import xml.etree.ElementTree as ET

MODELS_DIR = r"PROJECT_DIR\models\calibrated_devita_nlsr"

# Subject scale factors relative to S1 baseline
SCALE_FACTORS = {
    'S1': {'humerus': np.array([1.0000, 1.0000, 1.0000]), 'forearm': np.array([1.0000, 1.0000, 1.0000])},
    'S2': {'humerus': np.array([1.0357, 1.0357, 1.0357]), 'forearm': np.array([1.0456, 1.0456, 1.0456])},
    'S3': {'humerus': np.array([1.0501, 1.0501, 1.0501]), 'forearm': np.array([1.0984, 1.0984, 1.0984])},
    'S4': {'humerus': np.array([1.0003, 1.0003, 1.0003]), 'forearm': np.array([1.0812, 1.0812, 1.0812])},
    'S5': {'humerus': np.array([0.9598, 0.9598, 0.9598]), 'forearm': np.array([1.0605, 1.0605, 1.0605])},
    'S6': {'humerus': np.array([1.0375, 1.0375, 1.0375]), 'forearm': np.array([1.0386, 1.0386, 1.0386])},
}

# Baseline S1 coordinates (from anatomical axis optimization)
JOINT_CENTER = np.array([0.013144, -0.286273, -0.009595])

ORIGINS = {
    'cRCL': np.array([0.02047337, -0.28555095, 0.02198155]),
    'aRCL': np.array([0.02193452, -0.28555095, 0.02164239]),
    'pRCL': np.array([0.01901222, -0.28555095, 0.02232071]),
    'LUCL': np.array([0.01852516, -0.28555095, 0.02243376]),
}

INSERTIONS = {
    'cRCL': np.array([-0.00521398, -0.01434433, 0.03443540]), # on radius
    'aRCL': np.array([-0.00521398, -0.01434433, 0.03443540]), # on radius
    'pRCL': np.array([-0.00521398, -0.01434433, 0.03443540]), # on radius
    'LUCL': np.array([-0.02644393, -0.02926826, 0.01088904]), # on ulna
}

# Biophysical parameter definitions
PARAMS = {
    'LUCL': {'K': '49.200000', 'area': '13.35', 'eta': '4.9953', 'F_uts': '102.4', 'alpha': '1.01', 'beta': '0.080000'},
    'aRCL': {'K': '49.200000', 'area': '5.67',  'eta': '3.7570', 'F_uts': '36.67', 'alpha': '1.01', 'beta': '0.080000'},
    'cRCL': {'K': '49.200000', 'area': '5.67',  'eta': '3.7570', 'F_uts': '36.67', 'alpha': '1.01', 'beta': '0.080000'},
    'pRCL': {'K': '49.200000', 'area': '5.67',  'eta': '3.7570', 'F_uts': '36.67', 'alpha': '1.01', 'beta': '0.080000'},
    'tUCL': {'area': '2.23'},
    'ALCL': {'area': '4.25'},
}

def main():
    print("=================================================================")
    print("PHASE 1: UPDATING SUBJECT-SPECIFIC OPENSIM MODELS (S1 - S6)")
    print("=================================================================")
    
    for sid in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
        osim_path = os.path.join(MODELS_DIR, f"Elbow_model-{sid}.osim")
        tree = ET.parse(osim_path)
        root = tree.getroot()
        
        sh = SCALE_FACTORS[sid]['humerus']
        sf = SCALE_FACTORS[sid]['forearm']
        
        # Calculate subject-scaled attachment coordinates
        scaled_coords = {}
        for bname in ['cRCL', 'aRCL', 'pRCL', 'LUCL']:
            orig_scaled = JOINT_CENTER + (ORIGINS[bname] - JOINT_CENTER) * sh
            ins_scaled  = INSERTIONS[bname] * sf
            scaled_coords[bname] = {
                'origin': f"{orig_scaled[0]:.8f} {orig_scaled[1]:.8f} {orig_scaled[2]:.8f}",
                'insertion': f"{ins_scaled[0]:.8f} {ins_scaled[1]:.8f} {ins_scaled[2]:.8f}"
            }
            
        # Calculate subject-specific L0:
        # Scale L0 baseline with segment scaling:
        # S1 baseline lengths:
        # cRCL: L0 = 19.812 mm * scale
        # aRCL: L0 = 20.332 mm * scale (target eps_0 = 0.5% at extension)
        # pRCL: L0 = 20.321 mm * scale (target peak eps = 4.87% in flexion)
        # LUCL: L0 = 52.900 mm * scale (target eps_0 = 0.0% at extension)
        avg_scale = 0.5 * (np.mean(sh) + np.mean(sf))
        
        l0_sub = {
            'cRCL': 0.019812 * avg_scale,
            'aRCL': 0.020332 * avg_scale,
            'pRCL': 0.020321 * avg_scale,
            'LUCL': 0.052900 * avg_scale,
        }
        
        print(f"\n--- Participant {sid} (avg_scale={avg_scale:.4f}) ---")
        for bname in ['cRCL', 'aRCL', 'pRCL', 'LUCL']:
            print(f"  {bname}: L0 = {l0_sub[bname]*1000:.3f} mm")
            
        # Update <NonLinearLigament> elements
        for lig in root.findall('.//NonLinearLigament'):
            name = lig.get('name')
            if name in scaled_coords:
                # Update pathpoints
                pts = lig.findall('.//PathPoint')
                if len(pts) >= 2:
                    pts[0].find('location').text = scaled_coords[name]['origin']
                    pts[1].find('location').text = scaled_coords[name]['insertion']
                # Update L0
                lig.find('resting_length').text = f"{l0_sub[name]:.7f}"
                # Update params
                p = PARAMS[name]
                if lig.find('ka_a') is not None: lig.find('ka_a').text = p['K']
                if lig.find('areaa_a') is not None: lig.find('areaa_a').text = p['area']
                if lig.find('eta_a') is not None: lig.find('eta_a').text = p['eta']
                if lig.find('pcsa_force') is not None: lig.find('pcsa_force').text = p['F_uts']
                if lig.find('alphaa_a') is not None: lig.find('alphaa_a').text = p['alpha']
                if lig.find('beta_a') is not None: lig.find('beta_a').text = p['beta']
                
            elif name == 'tUCL':
                if lig.find('areaa_a') is not None: lig.find('areaa_a').text = PARAMS['tUCL']['area']
            elif name == 'ALCL':
                if lig.find('areaa_a') is not None: lig.find('areaa_a').text = PARAMS['ALCL']['area']
                
        # Update <Millard2012EquilibriumMuscle> mirror elements
        for mus in root.findall('.//Millard2012EquilibriumMuscle'):
            mname = mus.get('name')
            bname = mname.replace('_m', '')
            if bname in scaled_coords:
                pts = mus.findall('.//PathPoint')
                if len(pts) >= 2:
                    pts[0].find('location').text = scaled_coords[bname]['origin']
                    pts[1].find('location').text = scaled_coords[bname]['insertion']
                if mus.find('resting_length') is not None:
                    mus.find('resting_length').text = f"{l0_sub[bname]:.7f}"
                    
        # Write back updated .osim
        tree.write(osim_path, encoding='utf-8', xml_declaration=True)
        print(f"  [SUCCESS] Updated {osim_path}")

    print("\n=================================================================")
    print("PHASE 1 COMPLETE: All 6 subject models scaled and calibrated!")
    print("=================================================================")

if __name__ == '__main__':
    main()
