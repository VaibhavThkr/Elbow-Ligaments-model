import xml.etree.ElementTree as ET

t = ET.parse(r'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model\models\calibrated_devita_nlsr\Elbow_model-S1.osim').getroot()

ligs = ['cRCL', 'aRCL', 'pRCL', 'LUCL', 'tUCL', 'ALCL']
for lig in ligs:
    l = None
    for item in t.findall('.//NonLinearLigament'):
        if item.get('name') == lig:
            l = item
            break
    if l is not None:
        props = {c.tag: c.text.strip() for c in l if c.text and c.text.strip()}
        print(f"=== {lig} ===")
        for k in ['resting_length', 'pcsa_force', 'damping', 'stiffness_modulus', 'cross_sectional_area', 'exponential_scale_factor', 'exponential_coefficient']:
            print(f"  {k}: {props.get(k)}")
