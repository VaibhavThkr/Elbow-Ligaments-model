import xml.etree.ElementTree as ET

t = ET.parse(r'PROJECT_DIR\models\calibrated_devita_nlsr\Elbow_model-S1.osim').getroot()
for name in ['cRCL', 'aRCL', 'pRCL', 'LUCL', 'tUCL', 'ALCL']:
    for item in t.findall('.//NonLinearLigament'):
        if item.get('name') == name:
            print(f"=== {name} ===")
            for c in item:
                if c.tag != 'GeometryPath':
                    print(f"  <{c.tag}>{c.text}</{c.tag}>")
            break
