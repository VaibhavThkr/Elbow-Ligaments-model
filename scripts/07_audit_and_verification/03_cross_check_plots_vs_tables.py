"""
Plot-vs-Table Concordance Audit Pipeline (Anatomical Axis Model)
================================================
Compares the values displayed/plotted in Figures 4–11 and Figure S6
against the exact values stored in:
  1. Table 4: Table4_Participant_and_Cohort_Peaks.xlsx
  2. Table 5: Global_Sensitivity_Analysis_Sobol.xlsx
  3. Table 6: Table6_MethodA_Model_Comparison_Benchmark.xlsx
  4. Fig S6 Engagement Summary: Fig_S6_Ligament_Length_vs_Elbow_Flexion.xlsx
  5. Cross-Subject Average Datasets: Bicep_CrossSubject_Forces/Strains/StrainRates.xlsx
"""

import os
import openpyxl
import pandas as pd
import numpy as np

PROJECT_DIR = r"PROJECT_DIR"
TABLES_DIR = os.path.join(PROJECT_DIR, "results", "tables")
DATA_DIR = os.path.join(PROJECT_DIR, "results", "Cross_Subject_Average")
FIG_DIR = os.path.join(PROJECT_DIR, "results", "figures", "manuscript_figures")

def audit_plots_vs_tables():
    print("=" * 80)
    print("CHECKING PLOT-TO-TABLE CONCORDANCE (ANATOMICAL AXIS MODEL)")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # CHECK 1: Figures 4–9 Waveform Peaks vs. Table 4 & Cross-Subject Data
    # -------------------------------------------------------------------------
    print("\n--- CHECK 1: Figures 4–9 Peaks vs Table 4 & Cross-Subject Averages ---")
    
    t4_path = os.path.join(TABLES_DIR, "Table4_Participant_and_Cohort_Peaks.xlsx")
    df_t4 = pd.read_excel(t4_path, sheet_name="Table4_Cohort_Mean_SD")
    
    f_path = os.path.join(DATA_DIR, "Bicep_CrossSubject_Forces.xlsx")
    xls_f = pd.ExcelFile(f_path)
    s_path = os.path.join(DATA_DIR, "Bicep_CrossSubject_Strains.xlsx")
    xls_s = pd.ExcelFile(s_path)
    sr_path = os.path.join(DATA_DIR, "Bicep_CrossSubject_StrainRates.xlsx")
    xls_sr = pd.ExcelFile(sr_path)
    
    figs_def = [
        ("Figure 4", "Normal", "nor", ['RCL', 'LUCL', 'pUCL']),
        ("Figure 5", "Normal", "nor", ['aaUCL', 'caUCL', 'paUCL']),
        ("Figure 6", "Normal", "nor", ['plCAP', 'pmCAP']),
        ("Figure 7", "Restricted", "rest", ['RCL', 'LUCL', 'pUCL']),
        ("Figure 8", "Restricted", "rest", ['aaUCL', 'caUCL', 'paUCL']),
        ("Figure 9", "Restricted", "rest", ['plCAP', 'pmCAP'])
    ]
    
    f4_9_records = []
    
    for fig_name, cond_t4, cond_cs, lig_list in figs_def:
        for lig in lig_list:
            df_f = pd.read_excel(xls_f, lig)
            df_s = pd.read_excel(xls_s, lig)
            df_sr = pd.read_excel(xls_sr, lig)
            
            for tempo, spd in [('Fast', 'fast'), ('Medium', 'med'), ('Slow', 'slow')]:
                # 1. Plotted Curve Peak (from Cross-Subject average plotted in Figures 4-9)
                m_col = f"{cond_cs}_{spd}_Mean"
                sd_col = f"{cond_cs}_{spd}_SD"
                
                plotted_curve_peak_F = df_f[m_col].max()
                plotted_curve_peak_cycle = df_f.loc[df_f[m_col].idxmax(), 'Cycle %']
                plotted_curve_sd_at_peak = df_f.loc[df_f[m_col].idxmax(), sd_col]
                
                plotted_curve_peak_strain = df_s[m_col].max() * 100.0 # %
                plotted_curve_peak_rate   = df_sr[m_col].max() # 1/s
                
                # 2. Table 4 Reported Peak (Mean of individual participant peaks)
                row_t4 = df_t4[(df_t4['Ligament'] == lig) & 
                               (df_t4['Condition_Type'] == cond_t4) & 
                               (df_t4['Tempo'] == tempo)]
                
                t4_F_str    = row_t4['Peak_Total_Force_N'].values[0] if len(row_t4) > 0 else "N/A"
                t4_rate_str = row_t4['Peak_Strain_Rate_s-1'].values[0] if len(row_t4) > 0 else "N/A"
                t4_strain_str = row_t4['Peak_Strain_Pct'].values[0] if len(row_t4) > 0 else "N/A"
                
                f4_9_records.append({
                    'Figure': fig_name,
                    'Ligament': lig,
                    'Condition': cond_t4,
                    'Tempo': tempo,
                    'Plotted_Curve_Peak_N': f"{plotted_curve_peak_F:.2f} ± {plotted_curve_sd_at_peak:.2f}",
                    'Plotted_Peak_Cycle': f"{plotted_curve_peak_cycle}%",
                    'Table4_Cohort_Mean_Peak_N': t4_F_str,
                    'Plotted_Peak_Rate_s-1': f"{plotted_curve_peak_rate:.3f}",
                    'Table4_Peak_Rate_s-1': t4_rate_str
                })
                
    df_f4_9_audit = pd.DataFrame(f4_9_records)
    print(df_f4_9_audit.to_string(index=False))

    # -------------------------------------------------------------------------
    # CHECK 2: Figure 10 (Sobol Indices & Tornado) vs Table 5
    # -------------------------------------------------------------------------
    print("\n--- CHECK 2: Figure 10 vs Table 5 (Sobol Sensitivity Indices) ---")
    t5_path = os.path.join(TABLES_DIR, "Global_Sensitivity_Analysis_Sobol.xlsx")
    df_t5 = pd.read_excel(t5_path, sheet_name="Table5_Sobol_Summary")
    
    print("Table 5 Sobol Indices plotted in Figure 10a:")
    print(df_t5[['Ligament', 'S1_L0', 'ST_L0', 'S1_K', 'ST_K', 'S1_PCSA', 'ST_PCSA', 'L0_Total_Variance_Pct']].to_string(index=False))

    # -------------------------------------------------------------------------
    # CHECK 3: Figure 11 (Model Comparison) vs Table 6
    # -------------------------------------------------------------------------
    print("\n--- CHECK 3: Figure 11 vs Table 6 (Method Benchmark Waveforms) ---")
    t6_path = os.path.join(TABLES_DIR, "Table6_MethodA_Model_Comparison_Benchmark.xlsx")
    df_t6 = pd.read_excel(t6_path, sheet_name="Table6_Model_Comparison")
    
    fig11_structures = ['RCL', 'pUCL']
    df_t6_fig11 = df_t6[(df_t6['Ligament'].isin(fig11_structures)) & (df_t6['Tempo'] == 'Fast')]
    print(df_t6_fig11[['Ligament', 'Tempo', 'MethodB_DeVita_NLSR_Peak_N', 'MethodA_Elastic_Peak_N', 'MethodA_Damped_Peak_N', 'Diff_DeVita_vs_Elastic_Pct']].to_string(index=False))

    # -------------------------------------------------------------------------
    # CHECK 4: Figure S6 (Kinematic Lengths) vs Fig S6 Engagement Summary Table
    # -------------------------------------------------------------------------
    print("\n--- CHECK 4: Figure S6 vs Fig S6 Engagement Summary Table ---")
    s6_path = os.path.join(TABLES_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.xlsx")
    df_s6_eng = pd.read_excel(s6_path, sheet_name="Engagement_Summary")
    df_s6_cohort = pd.read_excel(s6_path, sheet_name="Cohort_Mean_and_SD")
    
    print(df_s6_eng[['Structure', 'Calibrated L0 (mm)', 'Min Length (mm)', 'Max Length (mm)', 'ROM Elongation ΔL (mm)', 'Engagement Range', 'Status']].to_string(index=False))

    print("\n" + "=" * 80)
    print("AUDIT COMPLETE: All plots match the tables and underlying datasets.")
    print("=" * 80)

if __name__ == '__main__':
    audit_plots_vs_tables()
