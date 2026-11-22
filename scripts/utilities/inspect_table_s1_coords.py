import xml.etree.ElementTree as ET

tree = ET.parse(r'd:\Documents\Elbow joint Project\03_Path_B_Anatomical_Axis_Model\models\calibrated_devita_nlsr\Elbow_model-S1.osim')
root = tree.getroot()

for lig in ['aRCL', 'cRCL', 'pRCL', 'LUCL', 'ALCL']:
    elem = root.find(f'.//NonLinearLigament[@name="{lig}"]')
    if elem is not None:
        print(f'=== {lig} ===')
        path = elem.find('.//GeometryPath')
        if path is not None:
            pts = path.findall('.//PathPoint')
            for pt in pts:
                name = pt.get('name')
                # Check socket or body
                body_elem = pt.find('body')
                body_str = body_elem.text if body_elem is not None else ''
                socket_elem = pt.find('socket_parent_frame')
                socket_str = socket_elem.text if socket_elem is not None else ''
                loc = pt.find('location').text if pt.find('location') is not None else ''
                frame = body_str or socket_str
                print(f'  Point {name}: frame="{frame}", location=[{loc}]')
            wraps = path.findall('.//PathWrap')
            for w in wraps:
                wo = w.find('wrap_object').text if w.find('wrap_object') is not None else ''
                method = w.find('method').text if w.find('method') is not None else ''
                print(f'  Wrap: object={wo}, method={method}')
