import os, glob, xml.etree.ElementTree as ET

for root, dirs, files in os.walk(r'BASE_DIR'):
    for f in files:
        if f.endswith('.xml'):
            fp = os.path.join(root, f)
            try:
                t = ET.parse(fp).getroot()
                scales = t.findall('.//Scale')
                if scales:
                    print(fp, len(scales))
                    for sc in scales:
                        seg = sc.find('segment_name').text if sc.find('segment_name') is not None else ''
                        factors = sc.find('scale_factors').text if sc.find('scale_factors') is not None else ''
                        print(f'   {seg}: {factors}')
            except Exception:
                pass
