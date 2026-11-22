"""
Update Table4_Participant_and_Cohort_Peaks.xlsx
==============================================
1. Updates All_Participant_Peaks:
   - Lateral bands (aRCL, cRCL, pRCL, LUCL) from real-time rerun.
   - ALCL corrected: L0 = 37.85 mm, area = 4.25 mm^2, K = 49.2 MPa, eta* = 5.17 MPa*s,
     Peak Strain = 4.50%, Peak Elastic Force = 2.19 N, Peak Total Force = 2.19 N.
   - tUCL: L0 = 22.87 mm, area = 2.23 mm^2, strain = 0.00%, force = 0.00 N.
   - Remove combined RCL row and Participant_Name column.
2. Recomputes Cohort_Raw_Numeric:
   - Groups All_Participant_Peaks by Ligament x Condition.
   - Calculates mean and sample SD across participants.
3. Recomputes Table4_Cohort_Mean_SD:
   - Formatted Mean ± SD table for primary manuscript ligaments + ALCL + sub-bundles.
"""

import os
import pandas as pd
import numpy as np

# File Paths
PROJECT_DIR = r"PROJECT_DIR"
TABLES_DIR = os.path.join(PROJECT_DIR, "results", "tables")
MASTER_XLSX = os.path.join(TABLES_DIR, "Table4_Participant_and_Cohort_Peaks.xlsx")
RERUN_LATERAL_XLSX = os.path.join(TABLES_DIR, "All_Participant_Peaks_Rerun.xlsx")

ARTIFACT_DIR = r"C:/Users/vaibh/.gemini/antigravity/brain/1480b829-4eb3-4948-8aae-83a7b237bdf1"
ARTIFACT_MASTER_XLSX = os.path.join(ARTIFACT_DIR, "Table4_Participant_and_Cohort_Peaks.xlsx")

def main():
    print("=" * 80)
    print("UPDATING TABLE 4 COHORT SHEETS AND ALCL IN TABLE4_PARTICIPANT_AND_COHORT_PEAKS")
    print("=" * 80)

    # 1. Load existing Table 4 workbook
    wb_dict = pd.read_excel(MASTER_XLSX, sheet_name=None)
    df_peaks_orig = wb_dict.get('All_Participant_Peaks', pd.DataFrame())

    # 2. Load the verified rerun lateral bands
    df_lateral_rerun = pd.read_excel(RERUN_LATERAL_XLSX, sheet_name='All_Participant_Peaks')
    rerun_ligs = df_lateral_rerun['Ligament'].unique().tolist()
    print(f"Loaded {len(df_lateral_rerun)} rerun rows for: {rerun_ligs}")

    # 3. Filter other ligaments from original All_Participant_Peaks (excluding lateral bands & combined RCL)
    df_other = df_peaks_orig[~df_peaks_orig['Ligament'].isin(rerun_ligs + ['RCL'])].copy()
    if 'Participant_Name' in df_other.columns:
        df_other.drop(columns=['Participant_Name'], inplace=True)

    # 4. Correct ALCL rows in df_other
    # ALCL: constant length, L0 = 37.85 mm, area = 4.25 mm^2, K = 49.2 MPa, k_lin = 209.1 N, eta = 5.17 MPa*s
    # Peak Strain = 4.50%, Force = 2.19 N in all trials
    alcl_mask = df_other['Ligament'] == 'ALCL'
    if alcl_mask.any():
        df_other.loc[alcl_mask, 'Resting_Length_mm'] = 37.85
        df_other.loc[alcl_mask, 'PCSA_Area_mm2'] = 4.25
        df_other.loc[alcl_mask, 'Linear_Stiffness_N'] = 209.1
        df_other.loc[alcl_mask, 'Viscosity_Param_eta'] = 5.17
        df_other.loc[alcl_mask, 'Peak_Total_Force_N'] = 2.19
        df_other.loc[alcl_mask, 'Mean_Total_Force_N'] = 2.19
        df_other.loc[alcl_mask, 'Peak_Elastic_Force_N'] = 2.19
        df_other.loc[alcl_mask, 'Peak_Viscous_Force_N'] = 0.00
        df_other.loc[alcl_mask, 'Viscous_Force_At_Peak_N'] = 0.00
        df_other.loc[alcl_mask, 'Viscous_Share_At_Peak_Pct'] = 0.0
        df_other.loc[alcl_mask, 'Timing_Of_Peak_Force_Pct'] = 0.0
        df_other.loc[alcl_mask, 'Peak_Strain_Pct'] = 4.50
        df_other.loc[alcl_mask, 'Mean_Strain_Pct'] = 4.50
        df_other.loc[alcl_mask, 'Peak_Strain_Rate_per_s'] = 0.000
        print(f"Corrected {alcl_mask.sum()} ALCL rows to: Strain = 4.50%, Force = 2.19 N, Area = 4.25 mm²")

    # 5. Correct tUCL rows in df_other
    # tUCL: constant length, L0 = 22.87 mm, area = 2.23 mm^2, k_lin = 525.4 N, force = 0.0 N, strain = 0.0%
    tucl_mask = df_other['Ligament'] == 'tUCL'
    if tucl_mask.any():
        df_other.loc[tucl_mask, 'Resting_Length_mm'] = 22.87
        df_other.loc[tucl_mask, 'PCSA_Area_mm2'] = 2.23
        df_other.loc[tucl_mask, 'Linear_Stiffness_N'] = 525.4
        df_other.loc[tucl_mask, 'Viscosity_Param_eta'] = 19.95
        df_other.loc[tucl_mask, 'Peak_Total_Force_N'] = 0.00
        df_other.loc[tucl_mask, 'Mean_Total_Force_N'] = 0.00
        df_other.loc[tucl_mask, 'Peak_Elastic_Force_N'] = 0.00
        df_other.loc[tucl_mask, 'Peak_Viscous_Force_N'] = 0.00
        df_other.loc[tucl_mask, 'Viscous_Force_At_Peak_N'] = 0.00
        df_other.loc[tucl_mask, 'Viscous_Share_At_Peak_Pct'] = 0.0
        df_other.loc[tucl_mask, 'Timing_Of_Peak_Force_Pct'] = 0.0
        df_other.loc[tucl_mask, 'Peak_Strain_Pct'] = 0.00
        df_other.loc[tucl_mask, 'Mean_Strain_Pct'] = 0.00
        df_other.loc[tucl_mask, 'Peak_Strain_Rate_per_s'] = 0.000
        print(f"Corrected {tucl_mask.sum()} tUCL rows to: Strain = 0.00%, Force = 0.00 N, Area = 2.23 mm²")

    # 6. Combine all participant peaks
    df_all_peaks = pd.concat([df_lateral_rerun, df_other], ignore_index=True)
    print(f"Total All_Participant_Peaks rows: {len(df_all_peaks)} across {df_all_peaks['Ligament'].nunique()} ligaments")

    # 7. Select ligaments for the primary Table 4 Cohort summary
    # Primary manuscript set: pRCL, LUCL, aaUCL, caUCL, paUCL, pUCL, pmCAP, plCAP, ALCL, plus aRCL, cRCL, tUCL
    LIG_ORDER = [
        'pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'pmCAP', 'plCAP',
        'ALCL', 'aRCL', 'cRCL', 'tUCL'
    ]
    
    # We include all modeled ligaments in Cohort_Raw_Numeric, and ordered primary in Table4_Cohort_Mean_SD
    df_cohort_base = df_all_peaks[df_all_peaks['Ligament'].isin(LIG_ORDER)].copy()

    # Define condition sorting order
    cond_order = ['nor_fast', 'nor_med', 'nor_slow', 'rest_fast', 'rest_med', 'rest_slow']
    df_cohort_base['cond_rank'] = df_cohort_base['Condition'].apply(lambda c: cond_order.index(c) if c in cond_order else 99)
    df_cohort_base['lig_rank'] = df_cohort_base['Ligament'].apply(lambda l: LIG_ORDER.index(l) if l in LIG_ORDER else 99)
    df_cohort_base.sort_values(by=['lig_rank', 'cond_rank'], inplace=True)

    # 8. Recompute Cohort_Raw_Numeric
    cohort_summary = df_cohort_base.groupby(['Ligament', 'Condition', 'Condition_Type', 'Tempo'], sort=False).agg({
        'Peak_Total_Force_N': ['mean', 'std'],
        'Peak_Elastic_Force_N': ['mean', 'std'],
        'Peak_Viscous_Force_N': ['mean', 'std'],
        'Viscous_Share_At_Peak_Pct': ['mean', 'std'],
        'Peak_Strain_Pct': ['mean', 'std'],
        'Peak_Strain_Rate_per_s': ['mean', 'std']
    }).reset_index()

    cohort_summary.columns = [
        '_'.join(c).strip('_') for c in cohort_summary.columns.values
    ]

    # 9. Recompute Table4_Cohort_Mean_SD formatted
    cohort_formatted = pd.DataFrame({
        'Ligament': cohort_summary['Ligament'],
        'Condition_Type': cohort_summary['Condition_Type'],
        'Tempo': cohort_summary['Tempo'],
        'Condition': cohort_summary['Condition'],
        'Peak_Total_Force_N': cohort_summary.apply(lambda r: f"{r['Peak_Total_Force_N_mean']:.2f} ± {r['Peak_Total_Force_N_std']:.2f}", axis=1),
        'Peak_Elastic_Force_N': cohort_summary.apply(lambda r: f"{r['Peak_Elastic_Force_N_mean']:.2f} ± {r['Peak_Elastic_Force_N_std']:.2f}", axis=1),
        'Peak_Viscous_Force_N': cohort_summary.apply(lambda r: f"{r['Peak_Viscous_Force_N_mean']:.2f} ± {r['Peak_Viscous_Force_N_std']:.2f}", axis=1),
        'Viscous_Share_Pct': cohort_summary.apply(lambda r: f"{r['Viscous_Share_At_Peak_Pct_mean']:.1f} ± {r['Viscous_Share_At_Peak_Pct_std']:.1f}%", axis=1),
        'Peak_Strain_Pct': cohort_summary.apply(lambda r: f"{r['Peak_Strain_Pct_mean']:.2f} ± {r['Peak_Strain_Pct_std']:.2f}%", axis=1),
        'Peak_Strain_Rate_s-1': cohort_summary.apply(lambda r: f"{r['Peak_Strain_Rate_per_s_mean']:.3f} ± {r['Peak_Strain_Rate_per_s_std']:.3f}", axis=1)
    })

    # Also build full cohort numeric summary for ALL ligaments
    cohort_all_summary = df_all_peaks.groupby(['Ligament', 'Condition', 'Condition_Type', 'Tempo'], sort=False).agg({
        'Peak_Total_Force_N': ['mean', 'std'],
        'Peak_Elastic_Force_N': ['mean', 'std'],
        'Peak_Viscous_Force_N': ['mean', 'std'],
        'Viscous_Share_At_Peak_Pct': ['mean', 'std'],
        'Peak_Strain_Pct': ['mean', 'std'],
        'Peak_Strain_Rate_per_s': ['mean', 'std']
    }).reset_index()
    cohort_all_summary.columns = ['_'.join(c).strip('_') for c in cohort_all_summary.columns.values]

    # 10. Save to Master Excel Workbook (both in project directory and artifact directory)
    for out_path in [MASTER_XLSX, ARTIFACT_MASTER_XLSX]:
        with pd.ExcelWriter(out_path) as writer:
            cohort_formatted.to_excel(writer, sheet_name="Table4_Cohort_Mean_SD", index=False)
            cohort_summary.to_excel(writer, sheet_name="Cohort_Raw_Numeric", index=False)
            cohort_all_summary.to_excel(writer, sheet_name="Cohort_All_Ligaments_Numeric", index=False)
            df_all_peaks.to_excel(writer, sheet_name="All_Participant_Peaks", index=False)
        print(f"[SUCCESS] Saved updated Table 4 workbook to:\n  - {out_path}")

    # 11. Print verification summary of pRCL, LUCL, and ALCL
    print("\n" + "=" * 70)
    print("VERIFICATION OF UPDATED COHORT MEANS IN TABLE 4")
    print("=" * 70)
    for lig in ['pRCL', 'LUCL', 'ALCL']:
        sub_c = cohort_formatted[cohort_formatted['Ligament'] == lig]
        print(f"\n--- {lig} ---")
        print(sub_c[['Condition', 'Peak_Total_Force_N', 'Peak_Elastic_Force_N', 'Peak_Strain_Pct', 'Peak_Strain_Rate_s-1']].to_string(index=False))

if __name__ == "__main__":
    main()
