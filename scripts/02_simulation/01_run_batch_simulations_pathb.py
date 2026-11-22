"""
Phase 3: Parallel Batch Forward Dynamic Simulations (Path B Anatomical Axis Model)
==================================================================================
Executes OpenSim AnalyzeTool (ForceReporter + MuscleAnalysis) across all 35 motion
files for all 6 subjects using the Path B calibrated models (Option A 3-fiber RCL fan)
and the compiled NLSR plugin, utilizing multi-threaded parallel execution across CPU cores.

Models:
  d:\\Documents\\Elbow joint Project\\03_Path_B_Anatomical_Axis_Model\\models\\calibrated_devita_nlsr\\Elbow_model-S{1..6}.osim

Outputs:
  d:\\Documents\\Elbow joint Project\\03_Path_B_Anatomical_Axis_Model\\simulations\\<Subject>\\Analysis Results\\<Trial>\\
    - <Trial>_ForceReporter_forces.sto
    - <Trial>_MuscleAnalysis_Length.sto
"""

import os
import glob
import subprocess
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

BASE_DIR = r"d:\Documents\Elbow joint Project"
PATHB_DIR = os.path.join(BASE_DIR, "03_Path_B_Anatomical_Axis_Model")
MODEL_DIR = os.path.join(PATHB_DIR, "models", "calibrated_devita_nlsr")
SIM_DIR = os.path.join(PATHB_DIR, "simulations")

OPENSIM_CMD = r"D:\Research_Workspace\Software_and_Tools\OpenSim 4.0\bin\opensim-cmd.exe"
PLUGIN_DLL  = r"D:\Research_Workspace\Software_and_Tools\OpenSim 4.0\plugins\NLSR_Plugin.dll"

# Subject mapping: (ID, Subject Name, Folder in Simulation and Analysis)
SUBJECTS = [
    ("S1", "Gokul",  "1_S1 Bicep"),
    ("S2", "Sanyam", "2_S2 Bicep"),
    ("S3", "Sunil",  "3_S3 Bicep"),
    ("S4", "Piyush", "4_S4 Bicep"),
    ("S5", "Nitin",  "5_S5 Bicep"),
    ("S6", "Rahul",  "6_S6 Bicep")
]

MAX_WORKERS = 6
print_lock = threading.Lock()

def get_mot_time_range(mot_path):
    with open(mot_path, "r") as f:
        lines = f.readlines()
    header_end = 0
    for i, line in enumerate(lines):
        if line.strip() == "endheader":
            header_end = i + 1
            break
    data_idx = header_end
    while data_idx < len(lines):
        parts = lines[data_idx].split()
        if parts:
            try:
                t0 = float(parts[0])
                break
            except ValueError:
                pass
        data_idx += 1
    last_row = lines[-1].split()
    return float(lines[data_idx].split()[0]), float(last_row[0])

def generate_analyze_xml(model_file, mot_file, out_dir, trial_name):
    start_t, end_t = get_mot_time_range(mot_file)
    xml_content = f"""<?xml version="1.0" encoding="UTF-8" ?>
<OpenSimDocument Version="40000">
	<AnalyzeTool name="{trial_name}">
		<model_file>{model_file}</model_file>
		<replace_force_set>false</replace_force_set>
		<force_set_files />
		<results_directory>{out_dir}</results_directory>
		<output_precision>8</output_precision>
		<initial_time>{start_t}</initial_time>
		<final_time>{end_t}</final_time>
		<AnalysisSet name="Analyses">
			<objects>
				<ForceReporter name="ForceReporter">
					<on>true</on>
					<start_time>{start_t}</start_time>
					<end_time>{end_t}</end_time>
					<step_interval>1</step_interval>
					<in_degrees>true</in_degrees>
					<include_forces>true</include_forces>
				</ForceReporter>
				<MuscleAnalysis name="MuscleAnalysis">
					<on>true</on>
					<start_time>{start_t}</start_time>
					<end_time>{end_t}</end_time>
					<step_interval>1</step_interval>
					<in_degrees>true</in_degrees>
					<muscle_list>all</muscle_list>
					<moment_arm_coordinate_list>none</moment_arm_coordinate_list>
					<compute_moments>false</compute_moments>
				</MuscleAnalysis>
			</objects>
		</AnalysisSet>
		<coordinates_file>{mot_file}</coordinates_file>
		<lowpass_cutoff_frequency_for_coordinates>6</lowpass_cutoff_frequency_for_coordinates>
	</AnalyzeTool>
</OpenSimDocument>"""
    setup_path = os.path.join(out_dir, "setup_analyze.xml")
    with open(setup_path, "w") as f:
        f.write(xml_content)
    return setup_path

def run_single_simulation(task_info):
    sid, name, folder_name, mot_file, idx, total = task_info
    model_file = os.path.join(MODEL_DIR, f"Elbow_model-{sid}.osim")
    trial_name = os.path.splitext(os.path.basename(mot_file))[0]
    out_dir = os.path.join(SIM_DIR, folder_name, "Analysis Results", trial_name)
    os.makedirs(out_dir, exist_ok=True)

    force_sto = os.path.join(out_dir, f"{trial_name}_ForceReporter_forces.sto")
    length_sto = os.path.join(out_dir, f"{trial_name}_MuscleAnalysis_Length.sto")

    # Check if already completed and valid (disabled when FORCE_RERUN is True)
    FORCE_RERUN = True
    if not FORCE_RERUN and os.path.exists(force_sto) and os.path.exists(length_sto):
        if os.path.getsize(force_sto) > 50000 and os.path.getsize(length_sto) > 50000:
            with print_lock:
                print(f"  [{idx:02d}/{total}] {sid} ({name}) -> {trial_name} [ALREADY EXISTS, SKIPPING]")
            return {
                "Subject_ID": sid,
                "Participant": name,
                "Trial": trial_name,
                "Duration_s": 0.0,
                "Status": "PASSED (CACHED)",
                "Force_STO_Exists": True,
                "Length_STO_Exists": True
            }

    with print_lock:
        print(f"  [{idx:02d}/{total}] STARTING {sid} ({name}) -> {trial_name} ...")

    setup_xml = generate_analyze_xml(model_file, mot_file, out_dir, trial_name)
    cmd = [OPENSIM_CMD, "-L", PLUGIN_DLL, "run-tool", setup_xml]

    t_start = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    duration = time.time() - t_start

    passed = (res.returncode == 0 and os.path.exists(force_sto) and os.path.exists(length_sto))

    with print_lock:
        if passed:
            print(f"  [{idx:02d}/{total}] FINISHED {sid} ({name}) -> {trial_name} [DONE ({duration:.1f}s)]")
        else:
            print(f"  [{idx:02d}/{total}] FAILED   {sid} ({name}) -> {trial_name} [FAILED ({duration:.1f}s)]")
            if res.stderr:
                print(f"    ERROR [{trial_name}]: {res.stderr.strip()[:200]}")

    return {
        "Subject_ID": sid,
        "Participant": name,
        "Trial": trial_name,
        "Duration_s": round(duration, 2),
        "Status": "PASSED" if passed else "FAILED",
        "Force_STO_Exists": os.path.exists(force_sto),
        "Length_STO_Exists": os.path.exists(length_sto)
    }

def main():
    print("=" * 75)
    print(f"STARTING PATH B PARALLEL FORWARD DYNAMIC SIMULATIONS ({MAX_WORKERS} Workers)")
    print("=" * 75)

    tasks = []
    task_idx = 0
    for sid, name, folder_name in SUBJECTS:
        mot_dir = os.path.join(BASE_DIR, "Simulation and Analysis", folder_name, "Motion files")
        mot_files = sorted(glob.glob(os.path.join(mot_dir, "*.mot")))
        for mf in mot_files:
            task_idx += 1
            tasks.append((sid, name, folder_name, mf, task_idx, 35))

    print(f"Total tasks gathered: {len(tasks)}")
    t0_overall = time.time()

    results_summary = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(run_single_simulation, t) for t in tasks]
        for f in as_completed(futures):
            results_summary.append(f.result())

    elapsed_total = time.time() - t0_overall
    success_count = sum(1 for r in results_summary if "PASSED" in r["Status"])

    print("\n" + "=" * 75)
    print(f"PATH B PARALLEL SIMULATION COMPLETE: {success_count}/{len(tasks)} successful.")
    print(f"Total wall time elapsed: {elapsed_total:.1f}s ({elapsed_total/60:.2f} min).")
    print(f"Outputs stored in: {SIM_DIR}")
    print("=" * 75)

    df_audit = pd.DataFrame(results_summary)
    df_audit = df_audit.sort_values(by=["Subject_ID", "Trial"])
    audit_path = os.path.join(PATHB_DIR, "results", "Phase3_Simulation_Batch_Audit.xlsx")
    df_audit.to_excel(audit_path, index=False)
    print(f"Audit log saved to: {audit_path}")

if __name__ == "__main__":
    main()
