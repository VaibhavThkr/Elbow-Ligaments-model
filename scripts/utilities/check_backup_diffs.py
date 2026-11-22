import xml.etree.ElementTree as ET

t1 = ET.parse(r'BASE_DIR\02_Earlier_Versions_Archive\legacy_models\Elbow_model-Gokul.osim.backup').getroot()
t2 = ET.parse(r'BASE_DIR\02_Earlier_Versions_Archive\legacy_models\Elbow_model-Sanyam.osim.backup').getroot()

diffs = []
for e1, e2 in zip(t1.iter(), t2.iter()):
    t1_txt = (e1.text or '').strip()
    t2_txt = (e2.text or '').strip()
    if e1.tag != e2.tag or e1.attrib != e2.attrib or t1_txt != t2_txt:
        diffs.append((e1.tag, e1.attrib, t1_txt, t2_txt))

print(f'Total diffs: {len(diffs)}')
for d in diffs[:20]:
    print(d)
