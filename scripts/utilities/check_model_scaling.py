import xml.etree.ElementTree as ET

for s in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
    p = rf'd:\Documents\Elbow joint Project\Simulation and Analysis\Opensim Elbow model\Elbow_model-{s}.osim'
    t = ET.parse(p).getroot()
    for m in t.findall('.//Millard2012EquilibriumMuscle'):
        if m.get('name') == 'biceplongus_r':
            pts = [(pt.get('name'), pt.find('location').text.strip()) for pt in m.findall('.//PathPoint')][:2]
            print(s, pts)
