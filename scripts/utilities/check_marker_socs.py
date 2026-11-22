import xml.etree.ElementTree as ET

for s in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
    p = rf'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model-{s}.osim'
    t = ET.parse(p).getroot()
    print(f'=== {s} ===')
    for m in t.findall('.//Marker'):
        name = m.get('name')
        if name in ['r_elb_l', 'r_elb_m', 'r_wri_l', 'r_wri_m']:
            loc = m.find('location').text.strip() if m.find('location') is not None else ''
            soc = m.find('socket_parent_frame').text.strip() if m.find('socket_parent_frame') is not None else ''
            print(f'  {name} ({soc}): {loc}')
