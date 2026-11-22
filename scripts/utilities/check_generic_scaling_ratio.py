import xml.etree.ElementTree as ET
import numpy as np

t_gen = ET.parse(r'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model (2).osim').getroot()
t_s1 = ET.parse(r'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model-S1.osim').getroot()

lig_names = ['aaUCL', 'caUCL', 'paUCL', 'pUCL', 'aAL', 'pAL', 'tUCL', 'amCAP', 'alCAP', 'pmCAP', 'plCAP']

for lig in lig_names:
    l_gen = None
    for item in t_gen.findall('.//NonLinearLigament'):
        if item.get('name') == lig: l_gen = item; break
    l_s1 = None
    for item in t_s1.findall('.//NonLinearLigament'):
        if item.get('name') == lig: l_s1 = item; break
    
    if l_gen is not None and l_s1 is not None:
        pts_g = l_gen.findall('.//PathPoint')
        pts_s = l_s1.findall('.//PathPoint')
        for pg, ps in zip(pts_g, pts_s):
            body = ps.find('socket_parent_frame').text if ps.find('socket_parent_frame') is not None else (ps.find('body').text if ps.find('body') is not None else '')
            cg = np.array([float(x) for x in pg.find('location').text.split()])
            cs = np.array([float(x) for x in ps.find('location').text.split()])
            ratio = cs / cg
            diff = cs - cg
            print(f'{lig} {ps.get("name")} ({body}): cg={cg}, cs={cs}, diff={diff}, ratio={ratio}')
