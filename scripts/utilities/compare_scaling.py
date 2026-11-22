import xml.etree.ElementTree as ET

for s in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
    p = rf'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model-{s}.osim'
    t = ET.parse(p).getroot()
    l = None
    for item in t.findall('.//NonLinearLigament'):
        if item.get('name') == 'aaUCL':
            l = item
            break
    pts = [(pt.get('name'), pt.find('location').text.strip()) for pt in l.findall('.//PathPoint')] if l is not None else []
    print(f'{s} aaUCL: {pts[0]}')
