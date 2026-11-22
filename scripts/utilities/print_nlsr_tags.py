import xml.etree.ElementTree as ET

t = ET.parse(r'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model\models\calibrated_devita_nlsr\Elbow_model-S1.osim').getroot()
for item in t.findall('.//NonLinearLigament'):
    if item.get('name') == 'LUCL':
        for c in item:
            print(f"<{c.tag}>{c.text}</{c.tag}>")
        break
