import xml.etree.ElementTree as ET
import numpy as np

def p3d(s):
    return np.array([float(x) for x in s.split()])

print("Participant Segment Lengths from Markers:")
for s in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
    p = rf'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model-{s}.osim'
    t = ET.parse(p).getroot()
    
    m_dict = {}
    for m in t.findall('.//Marker'):
        name = m.get('name')
        loc = m.find('location').text.strip() if m.find('location') is not None else '0 0 0'
        m_dict[name] = p3d(loc)
    
    # Humerus length: r_should to r_elb_m (or joint)
    # Humerus width: r_elb_l to r_elb_m
    if 'r_elb_l' in m_dict and 'r_elb_m' in m_dict:
        epicondyle_width = np.linalg.norm(m_dict['r_elb_l'] - m_dict['r_elb_m'])
        print(f'{s}: Epicondylar width = {epicondyle_width*1000:.2f} mm')
    
    if 'r_wri_l' in m_dict and 'r_wri_m' in m_dict:
        wrist_width = np.linalg.norm(m_dict['r_wri_l'] - m_dict['r_wri_m'])
        print(f'{s}: Wrist width = {wrist_width*1000:.2f} mm')
