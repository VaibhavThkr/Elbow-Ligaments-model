import xml.etree.ElementTree as ET

for s in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
    p = rf'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model\models\calibrated_devita_nlsr\Elbow_model-{s}.osim'
    t = ET.parse(p).getroot()
    print(f'=== {s} ===')
    for lig in ['cRCL', 'aRCL', 'pRCL', 'LUCL']:
        l = None
        for item in t.findall('.//NonLinearLigament'):
            if item.get('name') == lig:
                l = item
                break
        pts = [(pt.get('name'), pt.find('location').text.strip()) for pt in l.findall('.//PathPoint')] if l is not None else []
        l0 = l.find('resting_length').text.strip() if l is not None else 'None'
        print(f'  {lig} (L0={l0}): {pts}')
