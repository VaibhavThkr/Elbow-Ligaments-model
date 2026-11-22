import re

for s in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
    fp = rf'BASE_DIR\Simulation and Analysis\Opensim Elbow model\Elbow_model-{s}.osim'
    with open(fp, 'r', encoding='utf-8') as f:
        text = f.read()
    
    print(f'=== {s} ===')
    for lig in ['aaUCL', 'RCL', 'LUCL', 'pUCL']:
        m = re.search(r'<NonLinearLigament name="' + lig + r'">.*?</NonLinearLigament>', text, re.DOTALL)
        if m:
            pts = re.findall(r'<PathPoint name="(.*?)">.*?<location>(.*?)</location>', m.group(0), re.DOTALL)
            print(f'  {lig}: {pts}')
