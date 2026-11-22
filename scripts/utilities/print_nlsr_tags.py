import xml.etree.ElementTree as ET

t = ET.parse(r'PROJECT_DIR\models\calibrated_devita_nlsr\Elbow_model-S1.osim').getroot()
for item in t.findall('.//NonLinearLigament'):
    if item.get('name') == 'LUCL':
        for c in item:
            print(f"<{c.tag}>{c.text}</{c.tag}>")
        break
