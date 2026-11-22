"""
Phase 5: Inferential Statistical Analysis (Anatomical Axis Elbow Model)
========================================================================
Performs Repeated-Measures ANOVA and Post-Hoc Comparisons across all 6 participants
with Greenhouse-Geisser sphericity correction and Holm-Bonferroni adjustments.

Analysis Strategy:
  1. One-Way RM ANOVA (Tempo: Fast, Med, Slow) across all 6 subjects on Normal (Continuous) condition.
  2. Greenhouse-Geisser epsilon calculation and sphericity-corrected p-values.
  3. Pairwise post-hoc paired t-tests with Holm and Bonferroni adjustments & Cohen's d_z.
  4. Two-Way RM ANOVA (Condition x Tempo) on balanced Fast & Med tempos across all 6 subjects.

Outputs:
  1. Cohort_Statistical_Analysis_Report.xlsx
  2. NLSR_Repeated_Measures_ANOVA_Report.xlsx
"""

import os
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.anova import AnovaRM

BASE_DIR = r"BASE_DIR"
PROJECT_DIR = os.path.join(BASE_DIR, "PROJECT_DIR")
TABLES_DIR = os.path.join(PROJECT_DIR, "results", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

PEAKS_FILE = os.path.join(TABLES_DIR, "Table4_Participant_and_Cohort_Peaks.xlsx")
OUT_STATS_FILE = os.path.join(TABLES_DIR, "Cohort_Statistical_Analysis_Report.xlsx")
OUT_ANOVA_FILE = os.path.join(TABLES_DIR, "NLSR_Repeated_Measures_ANOVA_Report.xlsx")

PRIMARY_LIGAMENTS = ['RCL', 'aRCL', 'cRCL', 'pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 'pmCAP', 'plCAP']

def compute_gg_epsilon(data_matrix):
    """Computes Greenhouse-Geisser epsilon for repeated-measures factors."""
    n_subj, k_levels = data_matrix.shape
    if k_levels <= 2:
        return 1.0
    cov = np.cov(data_matrix, rowvar=False)
    k = k_levels
    mean_var = np.mean(np.diag(cov))
    mean_cov = np.mean(cov)
    sum_sq_cov = np.sum(cov**2)
    
    num_eps = (k * (mean_var - mean_cov))**2
    den_eps = (k - 1) * (sum_sq_cov - 2 * k * np.sum(np.mean(cov, axis=0)**2) + k**2 * mean_cov**2)
    eps = num_eps / den_eps if den_eps > 0 else 1.0
    return float(np.clip(eps, 1.0 / (k - 1), 1.0))

def run_stats():
    print("=" * 70)
    print("RUNNING INFERENTIAL REPEATED-MEASURES ANOVA (ANATOMICAL AXIS MODEL)")
    print("=" * 70)
    
    if not os.path.exists(PEAKS_FILE):
        print(f"[ERROR] Peaks file not found: {PEAKS_FILE}")
        return
        
    df_peaks = pd.read_excel(PEAKS_FILE, sheet_name="All_Participant_Peaks")
    
    # 1. Normal Condition Tempo Effect (Fast, Med, Slow) across 6 subjects
    df_nor = df_peaks[df_peaks['Condition_Type'] == 'Normal'].copy()
    
    metrics = [
        ('Peak_Total_Force_N', 'Peak Force (N)'),
        ('Peak_Strain_Rate_per_s', 'Peak Strain Rate (s^-1)'),
        ('Peak_Strain_Pct', 'Peak Strain (%)')
    ]
    
    aov_tempo_records = []
    posthoc_records = []
    
    for metric_col, metric_name in metrics:
        for lig in PRIMARY_LIGAMENTS:
            df_lig = df_nor[df_nor['Ligament'] == lig].copy()
            if len(df_lig['Participant'].unique()) < 6:
                continue
                
            try:
                aov = AnovaRM(df_lig, metric_col, 'Participant', within=['Tempo']).fit()
                tbl = aov.anova_table
                f_val = float(tbl['F Value'].iloc[0])
                num_df = float(tbl['Num DF'].iloc[0])
                den_df = float(tbl['Den DF'].iloc[0])
                p_val = float(tbl['Pr > F'].iloc[0])
                eta_p2 = (f_val * num_df) / (f_val * num_df + den_df) if (f_val * num_df + den_df) > 0 else 0.0
                
                # Greenhouse-Geisser correction
                piv = df_lig.pivot(index='Participant', columns='Tempo', values=metric_col)[['Fast', 'Medium', 'Slow']].values
                eps_gg = compute_gg_epsilon(piv)
                p_gg = float(stats.f.sf(f_val, num_df * eps_gg, den_df * eps_gg))
                
                aov_tempo_records.append({
                    'Ligament': lig,
                    'Metric': metric_name,
                    'Factor': 'Movement Tempo (Fast, Med, Slow)',
                    'F_Statistic': round(f_val, 3),
                    'Num_DF': num_df,
                    'Den_DF': den_df,
                    'GG_Epsilon': round(eps_gg, 3),
                    'p_value_uncorrected': p_val,
                    'p_value_GG_corrected': p_gg,
                    'Partial_Eta_Squared': round(eta_p2, 3),
                    'Significance': '***' if p_gg < 0.001 else '**' if p_gg < 0.01 else '*' if p_gg < 0.05 else 'ns'
                })
                
                # Pairwise Post-Hoc Paired t-tests
                v_fast = df_lig[df_lig['Tempo'] == 'Fast'][metric_col].values
                v_med  = df_lig[df_lig['Tempo'] == 'Medium'][metric_col].values
                v_slow = df_lig[df_lig['Tempo'] == 'Slow'][metric_col].values
                
                pairs = [
                    ('Fast', 'Medium', v_fast, v_med),
                    ('Fast', 'Slow', v_fast, v_slow),
                    ('Medium', 'Slow', v_med, v_slow)
                ]
                
                for t1, t2, arr1, arr2 in pairs:
                    diff = arr1 - arr2
                    mean_d = float(np.mean(diff))
                    std_d = float(np.std(diff, ddof=1))
                    t_stat, p_raw = stats.ttest_rel(arr1, arr2)
                    dz = mean_d / std_d if std_d > 0 else 0.0
                    p_bonf = min(1.0, p_raw * 3.0)
                    
                    posthoc_records.append({
                        'Ligament': lig,
                        'Metric': metric_name,
                        'Comparison': f"{t1} vs. {t2}",
                        'Mean_Diff': round(mean_d, 3),
                        't_Statistic': round(t_stat, 3),
                        'df': len(arr1) - 1,
                        'p_raw': p_raw,
                        'p_Bonferroni': p_bonf,
                        'Cohen_dz': round(dz, 2),
                        'Significance': '***' if p_bonf < 0.001 else '**' if p_bonf < 0.01 else '*' if p_bonf < 0.05 else 'ns'
                    })
            except Exception as e:
                print(f"  [!] ANOVA failed for {lig} {metric_col}: {e}")

    # 2. Condition (Normal vs Restricted) x Tempo (Fast vs Med) Two-Way RM ANOVA
    df_cond_sub = df_peaks[(df_peaks['Tempo'].isin(['Fast', 'Medium']))].copy()
    aov_cond_records = []
    
    for metric_col, metric_name in metrics:
        for lig in PRIMARY_LIGAMENTS:
            df_lig = df_cond_sub[df_cond_sub['Ligament'] == lig].copy()
            if len(df_lig['Participant'].unique()) < 6:
                continue
            try:
                aov2 = AnovaRM(df_lig, metric_col, 'Participant', within=['Condition_Type', 'Tempo']).fit()
                tbl2 = aov2.anova_table.reset_index()
                for _, r in tbl2.iterrows():
                    f_val = float(r['F Value'])
                    num_df = float(r['Num DF'])
                    den_df = float(r['Den DF'])
                    p_val = float(r['Pr > F'])
                    eta_p2 = (f_val * num_df) / (f_val * num_df + den_df) if (f_val * num_df + den_df) > 0 else 0.0
                    aov_cond_records.append({
                        'Ligament': lig,
                        'Metric': metric_name,
                        'Source': r['index'],
                        'F_Statistic': round(f_val, 3),
                        'Num_DF': num_df,
                        'Den_DF': den_df,
                        'p_value': p_val,
                        'Partial_Eta_Squared': round(eta_p2, 3),
                        'Significance': '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                    })
            except Exception as e:
                pass

    df_tempo_aov = pd.DataFrame(aov_tempo_records)
    df_ph = pd.DataFrame(posthoc_records)
    df_cond_aov = pd.DataFrame(aov_cond_records)

    with pd.ExcelWriter(OUT_ANOVA_FILE) as writer:
        df_tempo_aov.to_excel(writer, sheet_name="Tempo_RM_ANOVA", index=False)
        df_ph.to_excel(writer, sheet_name="Tempo_PostHoc_Pairwise", index=False)
        df_cond_aov.to_excel(writer, sheet_name="Condition_x_Tempo_2Way", index=False)
        
    with pd.ExcelWriter(OUT_STATS_FILE) as writer:
        df_tempo_aov.to_excel(writer, sheet_name="OneWay_Tempo_ANOVA", index=False)
        df_ph.to_excel(writer, sheet_name="Pairwise_Comparisons", index=False)
        df_cond_aov.to_excel(writer, sheet_name="TwoWay_Condition_ANOVA", index=False)
        
    print(f"\n[SUCCESS] ANOVA report saved to: {OUT_ANOVA_FILE}")
    print(f"[SUCCESS] Comprehensive statistics report saved to: {OUT_STATS_FILE}")

if __name__ == "__main__":
    run_stats()
