import xml.etree.ElementTree as ET

for s in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
    p = rf'd:\Documents\Elbow joint Project\Simulation and Analysis\Opensim Elbow model\Elbow_model-{s}.osim'
    t = ET.parse(p).getroot()
    print(f'=== {s} ===')
    for m in t.findall('.//Marker'):
        body = m.find('body').text if m.find('body') is not None else ''
        name = m.get('name')
        loc = m.find('location').text.strip() if m.find('location') is not None else ''
        if any(b in body for b in ['humerus', 'radius', 'ulna']):
            print(f'  {name} ({body}): {loc}')
