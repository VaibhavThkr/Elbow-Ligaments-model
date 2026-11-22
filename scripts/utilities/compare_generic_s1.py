import xml.etree.ElementTree as ET

t_gen = ET.parse(r'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model (2).osim').getroot()
t_s1 = ET.parse(r'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model-S1.osim').getroot()

ligs = ['aaUCL', 'caUCL', 'paUCL', 'pUCL', 'RCL', 'LUCL']
for lig in ligs:
    l_gen = None
    for item in t_gen.findall('.//NonLinearLigament'):
        if item.get('name') == lig:
            l_gen = item
            break
    l_s1 = None
    for item in t_s1.findall('.//NonLinearLigament'):
        if item.get('name') == lig:
            l_s1 = item
            break
    pts_gen = [(p.get('name'), p.find('location').text.strip()) for p in l_gen.findall('.//PathPoint')] if l_gen is not None else 'None'
    pts_s1 = [(p.get('name'), p.find('location').text.strip()) for p in l_s1.findall('.//PathPoint')] if l_s1 is not None else 'None'
    print(f'{lig}:')
    print('  generic:', pts_gen)
    print('  S1:     ', pts_s1)
