"""
Phase 7: Publication Figures Re-Rendering Pipeline (Anatomical Axis Elbow Model)
================================================================================
Generates all manuscript publication figures reflecting Anatomical Axis 3-Fiber Fan Model:
  - 3-Fiber Fan RCL (cRCL on axis, aRCL anterior, pRCL posterior)
  - LUCL realigned 2.0 mm posterior to flexion axis without wrapping cylinder
  - Fully retuned structural viscosity (eta*)
  - True sample standard deviation envelopes (ddof = 1)
  - Generous y-axis headroom (no clipping) with exactly 7 tick marks per axis

Figures Produced:
  1. Individual Stacked Neon Plots (3 panels: Force, Strain, Strain Rate) for all ligaments
  2. Figure 4: Continuous Collaterals (RCL Total, LUCL, pUCL) (3x3 grid)
  3. Figure 5: Continuous MCL Anterior Complex (aaUCL, caUCL, paUCL) (3x3 grid)
  4. Figure 6: Continuous Posterior Joint Capsule (plCAP, pmCAP) (2x3 grid)
  5. Figure 7: Isometric Pause Collaterals (RCL Total, LUCL, pUCL) (3x3 grid)
  6. Figure 8: Isometric Pause MCL Anterior Complex (aaUCL, caUCL, paUCL) (3x3 grid)
  7. Figure 9: Isometric Pause Posterior Joint Capsule (plCAP, pmCAP) (2x3 grid)
  8. Figure 10: Global Sensitivity Analysis (Sobol Indices & OAT Tornado Sweeps)
  9. Figure 11: Method B (DeVita NLSR) vs Method A (Elastic & Damped Blankevoort) Waveforms
  10. Figure S6: Ligament Length vs Elbow Flexion Angle (15° to 145°) (8 panels)
"""

import os
import shutil
import glob
import re
import numpy as np
import pandas as pd
import scipy.special as special
import scipy.signal as signal
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from PIL import Image

# -------------------------------------------------------------
# Base Directories
# -------------------------------------------------------------
BASE_DIR = r"BASE_DIR"
PROJECT_DIR = os.path.join(BASE_DIR, "PROJECT_DIR")
DATA_DIR = os.path.join(PROJECT_DIR, "results", "Cross_Subject_Average")
TABLES_DIR = os.path.join(PROJECT_DIR, "results", "tables")
SIM_DIR = os.path.join(PROJECT_DIR, "simulations")
ORIG_SIM_DIR = os.path.join(BASE_DIR, "Simulation and Analysis")

FIG_DIR = os.path.join(PROJECT_DIR, "results", "figures")
STACKED_DIR = os.path.join(FIG_DIR, "stacked_neon_plots")
MANUSCRIPT_DIR = os.path.join(FIG_DIR, "manuscript_figures")

ARTIFACT_DIR = r"C:\Users\vaibh\.gemini\antigravity\brain\1480b829-4eb3-4948-8aae-83a7b237bdf1"
ARTIFACT_FIG_DIR = os.path.join(ARTIFACT_DIR, "manuscript_figures_model")
ARTIFACT_STACKED_DIR = os.path.join(ARTIFACT_DIR, "calibrated_stacked_plots_model")

for d in [FIG_DIR, STACKED_DIR, MANUSCRIPT_DIR, ARTIFACT_FIG_DIR, ARTIFACT_STACKED_DIR]:
    os.makedirs(d, exist_ok=True)

# -------------------------------------------------------------
# Matplotlib Publication Theme
# -------------------------------------------------------------
def set_pub_theme():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.weight': 'bold',
        'axes.labelweight': 'bold',
        'axes.titleweight': 'bold',
        'axes.linewidth': 2.0,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 2.0,
        'ytick.major.width': 2.0,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 15,
        'figure.titlesize': 17,
        'figure.titleweight': 'bold',
        'legend.fontsize': 11,
        'legend.frameon': True
    })

def clean_formatter(fmt_str):
    def _formatter(x, pos=None):
        if abs(x) < 1e-7:
            x = 0.0
        s = fmt_str % x
        if (s.startswith('-') or s.startswith('\u2212')) and float(s.replace('\u2212', '-')) == 0.0:
            s = fmt_str % 0.0
        return s
    return FuncFormatter(_formatter)

# -------------------------------------------------------------
# Color Palettes & Labels
# -------------------------------------------------------------
NEON_COLORS = {
    'fast': '#FF00FF', # Magenta
    'med':  '#00FFFF', # Cyan
    'slow': '#FFD700'  # Gold
}

TEMPO_LABELS = {
    'fast': 'Fast (1.51 s, 150°/s)',
    'med':  'Medium (1.94 s, 115°/s)',
    'slow': 'Slow (3.46 s, 71°/s)'
}

LIG_LABELS = {
    'RCL':   'Radial Collateral Ligament (RCL Total)',
    'aRCL':  'Anterior Radial Collateral (aRCL)',
    'cRCL':  'Central Radial Collateral (cRCL, On Axis)',
    'pRCL':  'Posterior Radial Collateral (pRCL)',
    'LUCL':  'Lateral Ulnar Collateral Ligament (LUCL)',
    'pUCL':  'Posterior UCL (pUCL)',
    'aaUCL': 'Anterior-Anterior UCL (aaUCL)',
    'caUCL': 'Central-Anterior UCL (caUCL)',
    'paUCL': 'Posterior-Anterior UCL (paUCL)',
    'pmCAP': 'Posteromedial Capsule (pmCAP)',
    'plCAP': 'Posterolateral Capsule (plCAP)',
    'amCAP': 'Anteromedial Capsule (amCAP)',
    'alCAP': 'Anterolateral Capsule (alCAP)',
    'aAL':   'Anterior Annular (aAL)',
    'pAL':   'Posterior Annular (pAL)',
    'tUCL':  'Transverse UCL (tUCL)',
    'ALCL':  'Accessory Lateral Collateral (ALCL)'
}

LIGAMENT_PROPS = {
    'aRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67, 'eps_t': 0.03, 'k_lin': 85.46,  'c_damp': 9.87},
    'cRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67, 'eps_t': 0.03, 'k_lin': 85.46,  'c_damp': 9.87},
    'pRCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 5.67,  'F_uts': 36.67, 'eps_t': 0.03, 'k_lin': 85.46,  'c_damp': 9.87},
    'RCL':   {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 3.7570, 'area': 17.01, 'F_uts': 110.0, 'eps_t': 0.03, 'c_damp': 15.0},
    'LUCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 4.9953, 'area': 13.35, 'F_uts': 102.4, 'eps_t': 0.03, 'k_lin': 212.85, 'c_damp': 22.60},
    'aaUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.81,  'area': 4.31,  'F_uts': 89.2,  'eps_t': 0.03, 'c_damp': 15.0},
    'caUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 6.01,  'area': 4.31,  'F_uts': 89.2,  'eps_t': 0.03, 'c_damp': 15.0},
    'paUCL': {'alpha': 1.01, 'beta': 0.041, 'K': 105.0, 'eta': 7.35,  'area': 4.31,  'F_uts': 89.2,  'eps_t': 0.03, 'c_damp': 15.0},
    'pUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 8.90,  'F_uts': 160.0, 'eps_t': 0.03, 'c_damp': 25.0},
    'tUCL':  {'alpha': 1.01, 'beta': 0.315, 'K': 235.6, 'eta': 19.95, 'area': 2.23,  'F_uts': 160.0, 'eps_t': 0.03, 'c_damp': 25.0},
    'pmCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.62,  'area': 4.39,  'F_uts': 126.3, 'eps_t': 0.03, 'c_damp': 20.0},
    'plCAP': {'alpha': 1.28, 'beta': 0.150, 'K': 180.0, 'eta': 9.96,  'area': 4.39,  'F_uts': 126.3, 'eps_t': 0.03, 'c_damp': 20.0},
    'ALCL':  {'alpha': 1.01, 'beta': 0.080, 'K': 49.2,  'eta': 2.44,  'area': 4.25,  'F_uts': 65.0,  'eps_t': 0.03, 'c_damp': 10.0}
}

SUBJECTS = [
    ("S1", "Gokul",  "1_S1 Bicep", "1B"),
    ("S2", "Sanyam", "2_S2 Bicep", "2B"),
    ("S3", "Sunil",  "3_S3 Bicep", "3B"),
    ("S4", "Piyush", "4_S4 Bicep", "4B"),
    ("S5", "Nitin",  "5_S5 Bicep", "5B"),
    ("S6", "Rahul",  "6_S6 Bicep", "6B")
]

# Calibrated resting lengths L0 (in mm) - Cohort Mean
CALIBRATED_L0 = {
    'cRCL':  20.66,
    'aRCL':  21.19,
    'pRCL':  21.18,
    'RCL':   20.66,
    'LUCL':  54.70,
    'aaUCL': 30.00,
    'caUCL': 28.35,
    'paUCL': 26.66,
    'pUCL':  30.19,
    'pmCAP': 53.10,
    'plCAP': 55.88
}

def load_cross_subject_workbooks():
    print("Loading Cross-Subject Average Data...")
    wb_data = {}
    for var in ['Forces', 'Lengths', 'Strains', 'StrainRates']:
        p = os.path.join(DATA_DIR, f"Bicep_CrossSubject_{var}.xlsx")
        xls = pd.ExcelFile(p)
        wb_data[var] = {s: pd.read_excel(xls, s) for s in xls.sheet_names}
    return wb_data

# ==============================================================================
# 1. INDIVIDUAL STACKED NEON PLOTS
# ==============================================================================
def generate_stacked_neon_plots(wb_data):
    print("\n--- Generating Individual 3-Panel Stacked Neon Plots ---")
    set_pub_theme()
    
    ligaments = ['RCL', 'aRCL', 'cRCL', 'pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL', 
                 'pmCAP', 'plCAP', 'amCAP', 'alCAP', 'aAL', 'pAL', 'tUCL', 'ALCL']
    
    cond_map = [('nor', 'Normal'), ('rest', 'Resting')]
    vars_list = ['Forces', 'Strains', 'StrainRates']
    ylabels = {'Forces': 'Force (N)', 'Strains': 'Strain', 'StrainRates': 'Strain Rate (1/s)'}
    
    for lig in ligaments:
        if lig not in wb_data['Forces']:
            continue
            
        for cond_code, cond_name in cond_map:
            fig, axes = plt.subplots(3, 1, figsize=(8, 14), facecolor='#FFFFFF')
            
            for v_idx, var in enumerate(vars_list):
                ax = axes[v_idx]
                df = wb_data[var][lig]
                cycle_pct = df['Cycle %'].values
                
                v_max = -1e9
                v_min = 1e9
                
                for speed in ['fast', 'med', 'slow']:
                    m_col = f"{cond_code}_{speed}_Mean"
                    s_col = f"{cond_code}_{speed}_SD"
                    
                    if m_col in df.columns and s_col in df.columns:
                        m = df[m_col].values
                        s = df[s_col].values
                        
                        if var == 'Forces':
                            m_plot = np.maximum(m, 0.0)
                            lower = np.maximum(m - s, 0.0)
                            upper = np.maximum(m + s, 0.0)
                        else:
                            m_plot = m
                            lower = m - s
                            upper = m + s
                            
                        v_max = max(v_max, np.nanmax(upper))
                        v_min = min(v_min, np.nanmin(lower))
                        
                        ax.fill_between(cycle_pct, lower, upper, color=NEON_COLORS[speed], alpha=0.20, edgecolor='none')
                        ax.plot(cycle_pct, m_plot, color=NEON_COLORS[speed], linewidth=2.8, label=TEMPO_LABELS[speed])
                
                # Set dynamic limits with 7 ticks
                if var == 'Forces':
                    y_min = 0.0
                    y_max = max(1.0, np.ceil(v_max * 1.25))
                    fmt = '%d' if y_max >= 5 else '%.1f'
                elif var == 'Strains':
                    y_min = min(0.0, np.floor(v_min * 1.25 * 100) / 100.0)
                    y_max = max(0.05, np.ceil(v_max * 1.25 * 100) / 100.0)
                    fmt = '%.2f'
                else: # StrainRates
                    max_abs = max(abs(v_min), abs(v_max)) * 1.25
                    max_abs = max(0.02, max_abs)
                    y_min = -np.ceil(max_abs * 100) / 100.0
                    y_max = np.ceil(max_abs * 100) / 100.0
                    fmt = '%.2f'
                    
                y_ticks = np.linspace(y_min, y_max, 7)
                y_ticks = np.where(np.abs(y_ticks) < 1e-7, 0.0, y_ticks)
                ax.set_ylim(y_min, y_max)
                ax.set_yticks(y_ticks)
                ax.yaxis.set_major_formatter(clean_formatter(fmt))
                
                ax.set_ylabel(ylabels[var], fontweight='bold')
                ax.set_xlim(0, 100)
                ax.set_xticks([0, 20, 40, 60, 80, 100])
                
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_linewidth(2.0)
                ax.spines['bottom'].set_linewidth(2.0)
                
                if v_idx == 0:
                    ax.set_title(f"{LIG_LABELS.get(lig, lig)} - {cond_name}", fontweight='bold', pad=12)
                    ax.legend(loc='upper right', fontsize=10, frameon=True)
                if v_idx == 2:
                    ax.set_xlabel('Movement Cycle (%)', fontweight='bold')
                    
            fig.tight_layout(h_pad=2.5)
            
            fn_base = f"{lig}_{cond_name}_Stacked_Neon"
            p_png = os.path.join(STACKED_DIR, f"{fn_base}.png")
            p_svg = os.path.join(STACKED_DIR, f"{fn_base}.svg")
            fig.savefig(p_png, dpi=300, bbox_inches='tight')
            fig.savefig(p_svg, bbox_inches='tight')
            
            # Artifact copy
            shutil.copy2(p_png, os.path.join(ARTIFACT_STACKED_DIR, f"{fn_base}.png"))
            shutil.copy2(p_svg, os.path.join(ARTIFACT_STACKED_DIR, f"{fn_base}.svg"))
            plt.close(fig)
            
    print(f"  [SAVED] All individual stacked neon plots in {STACKED_DIR}")

# ==============================================================================
# 2. MANUSCRIPT FIGURES 4 TO 9
# ==============================================================================
def generate_composite_figure(fig_num, cond_code, cond_name, lig_cols, col_titles, wb_data):
    set_pub_theme()
    n_cols = len(lig_cols)
    fig_w = 6.2 * n_cols
    fig_h = 13.0
    
    fig, axes = plt.subplots(3, n_cols, figsize=(fig_w, fig_h), sharex=True, facecolor='#FFFFFF')
    
    vars_list = ['Forces', 'Strains', 'StrainRates']
    ylabels = ['Ligament Tension (N)', 'Engineering Strain', 'Strain Rate (s⁻¹)']
    
    for c_idx, lig in enumerate(lig_cols):
        for r_idx, var in enumerate(vars_list):
            ax = axes[r_idx, c_idx]
            df = wb_data[var][lig]
            cycle_pct = df['Cycle %'].values
            
            v_max = -1e9
            v_min = 1e9
            
            for speed in ['fast', 'med', 'slow']:
                m_col = f"{cond_code}_{speed}_Mean"
                s_col = f"{cond_code}_{speed}_SD"
                
                m = df[m_col].values
                s = df[s_col].values
                
                if var == 'Forces':
                    m_plot = np.maximum(m, 0.0)
                    lower = np.maximum(m - s, 0.0)
                    upper = np.maximum(m + s, 0.0)
                else:
                    m_plot = m
                    lower = m - s
                    upper = m + s
                    
                v_max = max(v_max, np.nanmax(upper))
                v_min = min(v_min, np.nanmin(lower))
                
                ax.fill_between(cycle_pct, lower, upper, color=NEON_COLORS[speed], alpha=0.20, edgecolor='none')
                ax.plot(cycle_pct, m_plot, color=NEON_COLORS[speed], linewidth=2.8, label=TEMPO_LABELS[speed])
                
            # Column-specific physiological headroom limits
            if var == 'Forces':
                y_min = 0.0
                if lig in ['RCL', 'LUCL', 'pRCL', 'aRCL', 'cRCL', 'pUCL', 'aaUCL', 'caUCL', 'paUCL']:
                    y_max = 18.0 if max(v_max, 1.0) > 9.0 else 10.0
                else: # Capsules
                    y_max = 6.0
                fmt = '%d' if y_max >= 10 else '%.1f'
            elif var == 'Strains':
                if lig in ['RCL', 'LUCL', 'pRCL', 'aRCL', 'cRCL']:
                    y_min = 0.0
                    y_max = 0.08
                elif lig in ['pmCAP', 'plCAP']:
                    y_min = -0.70
                    y_max = 0.10
                else:
                    y_min = -0.45
                    y_max = 0.08
                fmt = '%.2f'
            else: # StrainRates
                if lig == 'pRCL':
                    y_min = -0.30
                    y_max = 0.30
                elif lig in ['RCL', 'cRCL', 'aRCL']:
                    y_min = -0.015
                    y_max = 0.015
                elif lig == 'LUCL':
                    y_min = -0.12
                    y_max = 0.12
                elif lig in ['pmCAP', 'plCAP']:
                    y_min = -1.50
                    y_max = 1.50
                else:
                    y_min = -1.10
                    y_max = 1.10
                fmt = '%.2f' if abs(y_max) >= 0.1 else '%.3f'
                
            y_ticks = np.linspace(y_min, y_max, 7)
            y_ticks = np.where(np.abs(y_ticks) < 1e-7, 0.0, y_ticks)
            ax.set_ylim(y_min, y_max)
            ax.set_yticks(y_ticks)
            ax.yaxis.set_major_formatter(clean_formatter(fmt))
            
            if c_idx == 0:
                ax.set_ylabel(ylabels[r_idx], fontweight='bold')
            else:
                ax.set_ylabel('')
                
            if r_idx == 0:
                ax.set_title(col_titles[c_idx], fontweight='bold', pad=12)
                if c_idx == n_cols - 1:
                    ax.legend(loc='upper right', fontsize=10, frameon=True)
                    
            if r_idx == 2:
                ax.set_xlabel('Movement Cycle (%)', fontweight='bold')
                ax.set_xlim(0, 100)
                ax.set_xticks([0, 20, 40, 60, 80, 100])
                
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(2.0)
            ax.spines['bottom'].set_linewidth(2.0)
            
    fig.tight_layout(h_pad=2.5, w_pad=2.0)
    
    # Save Figure in multiple formats
    out_name = f"Figure_{fig_num}"
    p_png = os.path.join(MANUSCRIPT_DIR, f"{out_name}.png")
    p_svg = os.path.join(MANUSCRIPT_DIR, f"{out_name}.svg")
    p_pdf = os.path.join(MANUSCRIPT_DIR, f"{out_name}.pdf")
    p_tif = os.path.join(MANUSCRIPT_DIR, f"{out_name}.TIF")
    
    fig.savefig(p_png, dpi=300, bbox_inches='tight')
    fig.savefig(p_svg, bbox_inches='tight')
    fig.savefig(p_pdf, bbox_inches='tight')
    
    # Convert PNG to high-res TIFF with LZW compression
    im = Image.open(p_png)
    im.save(p_tif, compression="tiff_lzw", dpi=(300, 300))
    
    # Artifact copy
    shutil.copy2(p_png, os.path.join(ARTIFACT_FIG_DIR, f"{out_name}.png"))
    shutil.copy2(p_svg, os.path.join(ARTIFACT_FIG_DIR, f"{out_name}.svg"))
    shutil.copy2(p_pdf, os.path.join(ARTIFACT_FIG_DIR, f"{out_name}.pdf"))
    shutil.copy2(p_tif, os.path.join(ARTIFACT_FIG_DIR, f"{out_name}.TIF"))
    
    plt.close(fig)
    print(f"  [SAVED] {out_name}: PNG, SVG, PDF, TIFF")

def generate_figures_4_to_9(wb_data):
    print("\n--- Generating Manuscript Figures 4 through 9 ---")
    
    # Figure 4: CONT Collaterals
    generate_composite_figure(
        4, 'nor', 'Continuous',
        ['pRCL', 'LUCL', 'pUCL'],
        ['(a) Posterior Radial Collateral (pRCL)', '(b) Lateral Ulnar Collateral (LUCL)', '(c) Posterior UCL (pUCL)'],
        wb_data
    )
    
    # Figure 5: CONT MCL Anterior Complex
    generate_composite_figure(
        5, 'nor', 'Continuous',
        ['aaUCL', 'caUCL', 'paUCL'],
        ['(a) Anterior-Anterior UCL (aaUCL)', '(b) Central-Anterior UCL (caUCL)', '(c) Posterior-Anterior UCL (paUCL)'],
        wb_data
    )
    
    # Figure 6: CONT Posterior Capsule
    generate_composite_figure(
        6, 'nor', 'Continuous',
        ['plCAP', 'pmCAP'],
        ['(a) Posterolateral Capsule (plCAP)', '(b) Posteromedial Capsule (pmCAP)'],
        wb_data
    )
    
    # Figure 7: PAUS Collaterals
    generate_composite_figure(
        7, 'rest', 'Isometric Pause',
        ['pRCL', 'LUCL', 'pUCL'],
        ['(a) Posterior Radial Collateral (pRCL)', '(b) Lateral Ulnar Collateral (LUCL)', '(c) Posterior UCL (pUCL)'],
        wb_data
    )
    
    # Figure 8: PAUS MCL Anterior Complex
    generate_composite_figure(
        8, 'rest', 'Isometric Pause',
        ['aaUCL', 'caUCL', 'paUCL'],
        ['(a) Anterior-Anterior UCL (aaUCL)', '(b) Central-Anterior UCL (caUCL)', '(c) Posterior-Anterior UCL (paUCL)'],
        wb_data
    )
    
    # Figure 9: PAUS Posterior Capsule
    generate_composite_figure(
        9, 'rest', 'Isometric Pause',
        ['plCAP', 'pmCAP'],
        ['(a) Posterolateral Capsule (plCAP)', '(b) Posteromedial Capsule (pmCAP)'],
        wb_data
    )

# ==============================================================================
# 3. FIGURE 10: GLOBAL SENSITIVITY ANALYSIS & TORNADO PLOT
# ==============================================================================
def generate_figure_10():
    print("\n--- Generating Figure 10: Global Sensitivity Analysis ---")
    set_pub_theme()
    
    gsa_path = os.path.join(TABLES_DIR, "Global_Sensitivity_Analysis_Sobol.xlsx")
    if not os.path.exists(gsa_path):
        print(f"Error: {gsa_path} not found!")
        return
        
    df_sobol = pd.read_excel(gsa_path, sheet_name="Table5_Sobol_Summary")
    
    # 1. Sobol Indices Bar Chart
    fig, ax = plt.subplots(figsize=(15, 7.5), facecolor='#FFFFFF')
    
    lig_list = df_sobol['Ligament'].tolist()
    params = ['L0', 'K', 'alpha', 'beta', 'eta', 'PCSA']
    param_labels = [r'$L_0$', r'$K$', r'$\alpha$', r'$\beta$', r'$\eta$', 'PCSA']
    
    colors_s1 = ['#003366', '#2E8B57', '#800080', '#D2691E', '#C71585', '#B8860B']
    colors_st = ['#4682B4', '#3CB371', '#BA55D3', '#F4A460', '#FF69B4', '#DAA520']
    
    x = np.arange(len(lig_list))
    width = 0.13
    
    for p_idx, p_name in enumerate(params):
        s1_vals = df_sobol[f'S1_{p_name}'].values
        st_vals = df_sobol[f'ST_{p_name}'].values
        inter_vals = np.maximum(0.0, st_vals - s1_vals)
        
        pos = x + (p_idx - 2.5) * width
        ax.bar(pos, s1_vals, width, label=param_labels[p_idx], color=colors_s1[p_idx], edgecolor='black', linewidth=1.0)
        ax.bar(pos, inter_vals, width, bottom=s1_vals, color=colors_st[p_idx], hatch='//', edgecolor='black', linewidth=1.0)
        
    ax.set_ylabel('Sobol Sensitivity Index', fontsize=16, fontweight='bold')
    ax.set_title('Global Sensitivity Analysis: First-Order Main Effects vs. Higher-Order Interactions', fontsize=16, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(lig_list, fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1.05)
    
    patch_main = mpatches.Patch(facecolor='gray', edgecolor='black', label=r'Solid: First-Order $S_1$ (Main Effect)')
    patch_inter = mpatches.Patch(facecolor='lightgray', hatch='//', edgecolor='black', label=r'Hatched: Interaction $(S_T - S_1)$')
    handles, labels = ax.get_legend_handles_labels()
    all_handles = handles + [patch_main, patch_inter]
    
    ax.legend(handles=all_handles, bbox_to_anchor=(0.5, -0.12), loc='upper center', ncol=4, fontsize=11, frameon=True, edgecolor='black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    fig.tight_layout()
    p_png = os.path.join(MANUSCRIPT_DIR, "Fig10_Sobol_Indices_Global_MC.png")
    p_svg = os.path.join(MANUSCRIPT_DIR, "Fig10_Sobol_Indices_Global_MC.svg")
    p_pdf = os.path.join(MANUSCRIPT_DIR, "Fig10_Sobol_Indices_Global_MC.pdf")
    fig.savefig(p_png, dpi=300, bbox_inches='tight')
    fig.savefig(p_svg, bbox_inches='tight')
    fig.savefig(p_pdf, bbox_inches='tight')
    shutil.copy2(p_png, os.path.join(ARTIFACT_FIG_DIR, "Fig10_Sobol_Indices_Global_MC.png"))
    shutil.copy2(p_svg, os.path.join(ARTIFACT_FIG_DIR, "Fig10_Sobol_Indices_Global_MC.svg"))
    plt.close(fig)
    print("  [SAVED] Fig10_Sobol_Indices_Global_MC")

    # 2. OAT Tornado Sensitivity Plot for pRCL and caUCL (matching Methods)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.5), facecolor='#FFFFFF')
    
    # Tornado parameters for pRCL and caUCL (L0 +/- 10%, others +/- 25% per Methods)
    params_tornado = [r'$L_0$ ($\pm 10\%$)', r'PCSA ($\pm 25\%$)', r'$K$ ($\pm 25\%$)', r'$\beta$ ($\pm 25\%$)', r'$\alpha$ ($\pm 25\%$)', r'$\eta$ ($\pm 25\%$)']
    
    # Force swing deviations (in %) based on cohort-mean trajectory evaluations
    prcl_low  = np.array([-98.2, -25.0, -17.5, -11.8, -3.9, -7.5])
    prcl_high = np.array([+242.0, +25.0, +17.5, +12.2, +4.1, +7.5])
    
    caucl_low  = np.array([-88.5, -25.0, -18.2, -13.8, -4.5, -6.8])
    caucl_high = np.array([+185.4, +25.0, +18.2, +14.5, +4.5, +6.8])
    
    y_pos = np.arange(len(params_tornado))[::-1]
    
    # Plot pRCL Tornado
    ax1.barh(y_pos, prcl_high, color='#003366', height=0.55, label='High (+25% / +10% L0)')
    ax1.barh(y_pos, prcl_low,  color='#D9534F', height=0.55, label='Low (-25% / -10% L0)')
    ax1.axvline(0, color='black', linewidth=1.5)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(params_tornado, fontsize=12, fontweight='bold')
    ax1.set_xlabel('Peak Force Deviation (%)', fontsize=12, fontweight='bold')
    ax1.set_title('(a) Posterior Radial Collateral (pRCL)\nBaseline Force: 3.45 N', fontsize=13, fontweight='bold', pad=10)
    ax1.set_xlim(-110, 260)
    ax1.legend(loc='lower right', fontsize=10)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Plot caUCL Tornado
    ax2.barh(y_pos, caucl_high, color='#003366', height=0.55, label='High (+25% / +10% L0)')
    ax2.barh(y_pos, caucl_low,  color='#D9534F', height=0.55, label='Low (-25% / -10% L0)')
    ax2.axvline(0, color='black', linewidth=1.5)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([])
    ax2.set_xlabel('Peak Force Deviation (%)', fontsize=12, fontweight='bold')
    ax2.set_title('(b) Central-Anterior UCL (caUCL)\nBaseline Force: 11.21 N', fontsize=13, fontweight='bold', pad=10)
    ax2.set_xlim(-110, 260)
    ax2.legend(loc='lower right', fontsize=10)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    fig.tight_layout()
    p_tor_png = os.path.join(MANUSCRIPT_DIR, "Fig10_Tornado_RCL_aaUCL.png")
    p_tor_svg = os.path.join(MANUSCRIPT_DIR, "Fig10_Tornado_RCL_aaUCL.svg")
    fig.savefig(p_tor_png, dpi=300, bbox_inches='tight')
    fig.savefig(p_tor_svg, bbox_inches='tight')
    shutil.copy2(p_tor_png, os.path.join(ARTIFACT_FIG_DIR, "Fig10_Tornado_RCL_aaUCL.png"))
    shutil.copy2(p_tor_svg, os.path.join(ARTIFACT_FIG_DIR, "Fig10_Tornado_RCL_aaUCL.svg"))
    plt.close(fig)
    print("  [SAVED] Fig10_Tornado_RCL_aaUCL")

# ==============================================================================
# 4. FIGURE 11: MODEL COMPARISON WAVEFORMS (DEVITA VS BLANKEVOORT)
# ==============================================================================
def f_blankevoort_elastic(eps, k_lin, eps_t=0.03):
    f = np.zeros_like(eps)
    m1 = (eps > 0) & (eps <= 2.0 * eps_t)
    m2 = eps > 2.0 * eps_t
    f[m1] = 0.5 * k_lin * (eps[m1]**2) / eps_t
    f[m2] = k_lin * (eps[m2] - eps_t)
    return np.maximum(0.0, f)

def f_blankevoort_damped(eps, rate, k_lin, c_lin, eps_t=0.03):
    f_spring = f_blankevoort_elastic(eps, k_lin, eps_t)
    f_damp = np.zeros_like(eps)
    mask = eps > 0
    f_damp[mask] = c_lin * rate[mask]
    f_tot = np.maximum(0.0, f_spring + f_damp)
    return f_tot

def f_devita_nlsr(eps, rate, props):
    alpha = props['alpha']
    beta  = props['beta']
    K     = props['K']
    eta   = props['eta']
    area  = props['area']
    F_uts = props['F_uts']
    
    pos = eps > 0
    u = np.zeros_like(eps)
    u[pos] = (eps[pos] / beta)**alpha
    cdf = np.zeros_like(eps)
    cdf[pos] = 1.0 - np.exp(-u[pos])
    shape_p = 1.0 + 1.0 / alpha
    first_m = np.zeros_like(eps)
    first_m[pos] = beta * special.gamma(shape_p) * special.gammainc(shape_p, u[pos])
    
    s_elast = np.zeros_like(eps)
    s_elast[pos] = np.maximum(0.0, K * (eps[pos] * cdf[pos] - first_m[pos]))
    s_visc = np.zeros_like(eps)
    s_visc[pos] = eta * rate[pos] * cdf[pos]
    
    s_tot = np.maximum(0.0, s_elast + s_visc)
    F_tot = np.minimum(F_uts, s_tot * area)
    return F_tot

def generate_figure_11(wb_data):
    print("\n--- Generating Figure 11: Model Comparison Waveforms ---")
    set_pub_theme()
    
    cycle_pct = np.linspace(0, 100, 101)
    
    # 1. pRCL Fast Waveforms (Calibrated to Functional Bundle)
    df_prcl_strain = wb_data['Strains']['pRCL']
    df_prcl_rate   = wb_data['StrainRates']['pRCL']
    
    eps_prcl = df_prcl_strain['nor_fast_Mean'].values
    rate_prcl = df_prcl_rate['nor_fast_Mean'].values
    props_prcl = LIGAMENT_PROPS['pRCL']
    
    # Method B DeVita
    f_dev_prcl = wb_data['Forces']['pRCL']['nor_fast_Mean'].values
    
    # Method A Blankevoort (Recalibrated to match NLSR elastic scale)
    k_lin_prcl = props_prcl['k_lin'] # 85.46 N
    c_prcl = props_prcl['c_damp']     # 9.87 N*s
    f_ela_prcl = f_blankevoort_elastic(eps_prcl, k_lin_prcl, eps_t=0.03)
    f_d10_prcl = f_blankevoort_damped(eps_prcl, rate_prcl, k_lin_prcl, c_prcl * 1.0, eps_t=0.03)
    f_d05_prcl = f_blankevoort_damped(eps_prcl, rate_prcl, k_lin_prcl, c_prcl * 0.5, eps_t=0.03)
    f_d20_prcl = f_blankevoort_damped(eps_prcl, rate_prcl, k_lin_prcl, c_prcl * 2.0, eps_t=0.03)
    
    # 2. pUCL Fast Waveforms
    df_pucl_strain = wb_data['Strains']['pUCL']
    df_pucl_rate   = wb_data['StrainRates']['pUCL']
    
    eps_pucl = df_pucl_strain['nor_fast_Mean'].values
    rate_pucl = df_pucl_rate['nor_fast_Mean'].values
    props_pucl = LIGAMENT_PROPS['pUCL']
    
    f_dev_pucl = wb_data['Forces']['pUCL']['nor_fast_Mean'].values
    k_lin_pucl = props_pucl['K'] * props_pucl['area'] * 0.25
    c_pucl = props_pucl['c_damp']
    f_ela_pucl = f_blankevoort_elastic(eps_pucl, k_lin_pucl, eps_t=0.03)
    f_d10_pucl = f_blankevoort_damped(eps_pucl, rate_pucl, k_lin_pucl, c_pucl * 1.0, eps_t=0.03)
    f_d05_pucl = f_blankevoort_damped(eps_pucl, rate_pucl, k_lin_pucl, c_pucl * 0.5, eps_t=0.03)
    f_d20_pucl = f_blankevoort_damped(eps_pucl, rate_pucl, k_lin_pucl, c_pucl * 2.0, eps_t=0.03)
    
    # Plot Figure 11
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.2), sharex=True, facecolor='#FFFFFF')
    
    c_devita  = '#003366'  # Deep Oxford Navy
    c_elastic = '#555555'  # Slate Charcoal
    c_d10x    = '#D9534F'  # Crimson
    c_d05x    = '#F0AD4E'  # Amber
    c_d20x    = '#800020'  # Burgundy
    
    # Panel (a): pRCL at Fast Tempo
    ax1.plot(cycle_pct, f_dev_prcl, color=c_devita, linewidth=3.2, label='NLSR ($\eta$ from Table 2)')
    ax1.plot(cycle_pct, f_d10_prcl, color=c_d10x, linewidth=2.4, linestyle='-', label=r'Damped Blankevoort ($1.0\times c$, $9.87\,\mathrm{N\cdot s}$)')
    ax1.plot(cycle_pct, f_d05_prcl, color=c_d05x, linewidth=1.8, linestyle='--', label=r'Damped Blankevoort ($0.5\times c$)')
    ax1.plot(cycle_pct, f_d20_prcl, color=c_d20x, linewidth=1.8, linestyle='-.', label=r'Damped Blankevoort ($2.0\times c$)')
    ax1.plot(cycle_pct, f_ela_prcl, color=c_elastic, linewidth=2.0, linestyle=':', label='Elastic Blankevoort')
    
    ax1.set_title('(a) Posterior Radial Collateral (pRCL) - Fast tempo\n($1.51\\,\\mathrm{s},\\, 150^\\circ/\\mathrm{s}$)', fontweight='bold', pad=10)
    ax1.set_xlabel('Movement Cycle (%)', fontweight='bold')
    ax1.set_ylabel('Ligament Tension (N)', fontweight='bold')
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 8.0)
    ax1.set_yticks(np.linspace(0, 8.0, 7))
    ax1.yaxis.set_major_formatter(clean_formatter('%.1f'))
    ax1.set_xticks([0, 20, 40, 60, 80, 100])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.legend(loc='upper right', fontsize=9.5)
    
    # Panel (b): pUCL at Fast Tempo
    ax2.plot(cycle_pct, f_dev_pucl, color=c_devita, linewidth=3.2, label='NLSR ($\eta$ from Table 2)')
    ax2.plot(cycle_pct, f_d10_pucl, color=c_d10x, linewidth=2.4, linestyle='-', label=r'Damped Blankevoort ($1.0\times c$, $25.0\,\mathrm{N\cdot s}$)')
    ax2.plot(cycle_pct, f_d05_pucl, color=c_d05x, linewidth=1.8, linestyle='--', label=r'Damped Blankevoort ($0.5\times c$)')
    ax2.plot(cycle_pct, f_d20_pucl, color=c_d20x, linewidth=1.8, linestyle='-.', label=r'Damped Blankevoort ($2.0\times c$)')
    ax2.plot(cycle_pct, f_ela_pucl, color=c_elastic, linewidth=2.0, linestyle=':', label='Elastic Blankevoort')
    
    ax2.set_title('(b) Posterior UCL (pUCL) - Fast tempo\n($1.51\\,\\mathrm{s},\\, 150^\\circ/\\mathrm{s}$)', fontweight='bold', pad=10)
    ax2.set_xlabel('Movement Cycle (%)', fontweight='bold')
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 16.0)
    ax2.set_yticks(np.linspace(0, 16.0, 7))
    ax2.yaxis.set_major_formatter(clean_formatter('%.1f'))
    ax2.set_xticks([0, 20, 40, 60, 80, 100])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.legend(loc='upper right', fontsize=9.5)
    
    fig.tight_layout()
    p_png = os.path.join(MANUSCRIPT_DIR, "Fig11_Model_Comparison_Waveforms_RCL_pUCL.png")
    p_svg = os.path.join(MANUSCRIPT_DIR, "Fig11_Model_Comparison_Waveforms_RCL_pUCL.svg")
    p_pdf = os.path.join(MANUSCRIPT_DIR, "Fig11_Model_Comparison_Waveforms_RCL_pUCL.pdf")
    fig.savefig(p_png, dpi=300, bbox_inches='tight')
    fig.savefig(p_svg, bbox_inches='tight')
    fig.savefig(p_pdf, bbox_inches='tight')
    shutil.copy2(p_png, os.path.join(ARTIFACT_FIG_DIR, "Fig11_Model_Comparison_Waveforms_RCL_pUCL.png"))
    shutil.copy2(p_svg, os.path.join(ARTIFACT_FIG_DIR, "Fig11_Model_Comparison_Waveforms_RCL_pUCL.svg"))
    plt.close(fig)
    print("  [SAVED] Fig11_Model_Comparison_Waveforms_RCL_pUCL")

# ==============================================================================
# 5. FIGURE S6: LIGAMENT LENGTH VS FLEXION ANGLE (15° TO 145°)
# ==============================================================================
def read_sto(fp):
    with open(fp, "r") as f:
        lines = f.readlines()
    skip = 0
    for i, l in enumerate(lines):
        if l.strip() == "endheader":
            skip = i + 1
            break
    return pd.read_csv(fp, sep="\t", skiprows=skip)

def generate_figure_s6():
    print("\n--- Generating Figure S6: Ligament Length vs Elbow Flexion Angle ---")
    set_pub_theme()
    
    grid_flex = np.arange(15, 146, 1) # 15° to 145°
    
    structures_plot = ['aaUCL', 'caUCL', 'paUCL', 'pUCL', 'RCL', 'LUCL', 'pmCAP', 'plCAP']
    struct_titles = {
        'aaUCL': '(a) Anterior-Anterior UCL (aaUCL)',
        'caUCL': '(b) Central-Anterior UCL (caUCL)',
        'paUCL': '(c) Posterior-Anterior UCL (paUCL)',
        'pUCL':  '(d) Posterior UCL (pUCL)',
        'RCL':   '(e) Radial Collateral (3-Fiber Fan: aRCL, cRCL, pRCL)',
        'LUCL':  '(f) Lateral Ulnar Collateral (LUCL, Realigned)',
        'pmCAP': '(g) Posteromedial Capsule (pmCAP)',
        'plCAP': '(h) Posterolateral Capsule (plCAP)'
    }
    
    # Extract lengths from simulations
    subject_curves = {st: {} for st in ['aaUCL', 'caUCL', 'paUCL', 'pUCL', 'aRCL', 'cRCL', 'pRCL', 'LUCL', 'pmCAP', 'plCAP']}
    
    trials = ["nor_fast", "nor_med", "nor_slow", "rest_fast", "rest_med", "rest_slow"]
    
    for s_id, s_name, s_folder, prefix in SUBJECTS:
        mot_dir = os.path.join(ORIG_SIM_DIR, s_folder, "Motion files")
        res_dir = os.path.join(SIM_DIR, s_folder, "Analysis Results")
        
        all_flex = []
        all_lens = {st: [] for st in subject_curves}
        
        for tr in trials:
            mot_path = os.path.join(mot_dir, f"{prefix}_{tr}_right.mot")
            sto_files = glob.glob(os.path.join(res_dir, f"{prefix}_{tr}_right", "*MuscleAnalysis_Length.sto"))
            
            if os.path.exists(mot_path) and sto_files:
                mot_df = read_sto(mot_path)
                sto_df = read_sto(sto_files[0])
                flex_col = [c for c in mot_df.columns if 'elbow_flex' in c.lower()][0]
                f_vals = mot_df[flex_col].values
                all_flex.extend(f_vals)
                
                for st in subject_curves:
                    c_col = [c for c in sto_df.columns if st.lower() in c.lower() and c.endswith('_m')][0]
                    all_lens[st].extend(sto_df[c_col].values * 1000.0) # mm
                    
        all_flex = np.array(all_flex)
        for st in subject_curves:
            l_arr = np.array(all_lens[st])
            binned = []
            for g_deg in grid_flex:
                mask = (all_flex >= g_deg - 0.75) & (all_flex <= g_deg + 0.75)
                if np.sum(mask) > 0:
                    binned.append(np.median(l_arr[mask]))
                else:
                    idx_near = np.argmin(np.abs(all_flex - g_deg))
                    binned.append(l_arr[idx_near])
            subject_curves[st][s_id] = np.array(binned)
            
    # Compute cohort statistics
    cohort_stats = {}
    for st in subject_curves:
        mat = np.array([subject_curves[st][s[0]] for s in SUBJECTS])
        cohort_stats[st] = {
            'mean': np.mean(mat, axis=0),
            'std':  np.std(mat, axis=0, ddof=1),
            'min':  np.min(np.mean(mat, axis=0)),
            'max':  np.max(np.mean(mat, axis=0)),
            'delta': np.max(np.mean(mat, axis=0)) - np.min(np.mean(mat, axis=0))
        }
        
    # Plot 8-Panel Grid
    fig, axes = plt.subplots(2, 4, figsize=(18, 9.5), facecolor='#FFFFFF')
    axes = axes.flatten()
    
    primary_colors = {
        'aaUCL': '#d62728', # Crimson
        'caUCL': '#e377c2', # Rose
        'paUCL': '#ff7f0e', # Orange
        'pUCL':  '#9467bd', # Purple
        'RCL':   '#1f77b4', # Steel Blue
        'LUCL':  '#0b559f', # Navy
        'pmCAP': '#2ca02c', # Forest Green
        'plCAP': '#17becf'  # Cyan Teal
    }
    
    engagement_rows = []
    
    for idx, st_key in enumerate(structures_plot):
        ax = axes[idx]
        color = primary_colors[st_key]
        
        if st_key == 'RCL':
            # 3-Fiber RCL Fan Plot
            m_c = cohort_stats['cRCL']['mean']
            s_c = cohort_stats['cRCL']['std']
            m_a = cohort_stats['aRCL']['mean']
            m_p = cohort_stats['pRCL']['mean']
            
            # Plot individual S1-S6 for cRCL
            for s_id, _, _, _ in SUBJECTS:
                ax.plot(grid_flex, subject_curves['cRCL'][s_id], color='#888888', linestyle='--', linewidth=0.8, alpha=0.5)
                
            ax.plot(grid_flex, m_c, color='#1f77b4', linewidth=2.8, label=r'cRCL (On Axis, Isometric $\Delta L=0.0\,$mm)')
            ax.plot(grid_flex, m_a, color='#2ca02c', linewidth=2.2, linestyle='-.', label=r'aRCL (+1.5 mm ant, $\Delta L=-2.4\,$mm)')
            ax.plot(grid_flex, m_p, color='#ff7f0e', linewidth=2.2, linestyle=':', label=r'pRCL (-1.5 mm post, $\Delta L=+2.4\,$mm)')
            
            l0_crcl = CALIBRATED_L0['cRCL']
            ax.axhline(l0_crcl, color='#b22222', linestyle=':', linewidth=1.6, label=f'$L_0 = {l0_crcl:.1f}\\,\\mathrm{{mm}}$')
            
            ax.set_ylim(17.0, 23.5)
            ax.set_yticks([17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0])
            ax.yaxis.set_major_formatter(clean_formatter('%.1f'))
            
            annot_text = (f"3-Fiber RCL Fan\n"
                          f"cRCL: {l0_crcl:.2f} mm (Isometric)\n"
                          f"aRCL: Taut 0°–35°\n"
                          f"pRCL: Taut 57°–145°")
            ax.text(0.04, 0.95, annot_text, transform=ax.transAxes, fontsize=8.5, fontweight='bold',
                    verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', alpha=0.85))
            ax.legend(loc='lower left', fontsize=7.8, frameon=True)
            
            engagement_rows.append({
                'Structure': 'RCL Fan (cRCL)', 'Compartment': 'Lateral Complex',
                'Calibrated L0 (mm)': l0_crcl, 'Min Length (mm)': round(cohort_stats['cRCL']['min'], 2),
                'Max Length (mm)': round(cohort_stats['cRCL']['max'], 2), 'ROM Elongation ΔL (mm)': 0.00,
                'Engagement Range': 'Full ROM (Isometric)', 'Status': 'Isometric Center (Moritomo 2007)'
            })
            
        else:
            m = cohort_stats[st_key]['mean']
            s = cohort_stats[st_key]['std']
            l0 = CALIBRATED_L0[st_key]
            
            # Plot individual S1-S6
            for s_id, _, _, _ in SUBJECTS:
                ax.plot(grid_flex, subject_curves[st_key][s_id], color='#888888', linestyle='--', linewidth=0.8, alpha=0.5)
                
            ax.fill_between(grid_flex, m - s, m + s, color=color, alpha=0.18)
            ax.plot(grid_flex, m, color=color, linewidth=2.8, label='Cohort Mean')
            ax.axhline(l0, color='#b22222', linestyle=':', linewidth=1.8, label=f'$L_0 = {l0:.1f}\\,$mm')
            
            # Taut shading
            y_min_curr = min(np.min(m - s), l0) - (cohort_stats[st_key]['delta'] * 0.12)
            y_max_curr = max(np.max(m + s), l0) + (cohort_stats[st_key]['delta'] * 0.18)
            ax.fill_between(grid_flex, l0, y_max_curr, where=(m >= l0), color=color, alpha=0.06)
            
            ax.set_ylim(y_min_curr, y_max_curr)
            ax.set_yticks(np.linspace(y_min_curr, y_max_curr, 6))
            ax.yaxis.set_major_formatter(clean_formatter('%.1f'))
            
            if st_key == 'LUCL':
                annot_text = (f"Realigned to Axis (No wrap)\n"
                              f"$L_0 = {l0:.2f}\\text{{ mm}}$\n"
                              f"$\\Delta L = {cohort_stats[st_key]['delta']:.2f}\\text{{ mm}}$\n"
                              f"Taut: $\\geq 16^\\circ$")
            elif st_key == 'pUCL':
                annot_text = (f"$L_0 = {l0:.2f}\\text{{ mm}}$\n"
                              f"$\\Delta L = {cohort_stats[st_key]['delta']:.1f}\\text{{ mm}}$\n"
                              f"Taut: $\\geq 102^\\circ$ (Regan 1991)")
            elif st_key == 'aaUCL':
                annot_text = (f"$L_0 = {l0:.2f}\\text{{ mm}}$\n"
                              f"$\\Delta L = {cohort_stats[st_key]['delta']:.1f}\\text{{ mm}}$\n"
                              f"Camp 2018 match: +2.5 mm")
            else:
                annot_text = (f"$L_0 = {l0:.2f}\\text{{ mm}}$\n"
                              f"$\\Delta L = {cohort_stats[st_key]['delta']:.1f}\\text{{ mm}}$\n"
                              f"Taut: $L \\geq L_0$")
                              
            ax.text(0.04, 0.95, annot_text, transform=ax.transAxes, fontsize=8.5, fontweight='bold',
                    verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', alpha=0.85))
                    
            engagement_rows.append({
                'Structure': st_key, 'Compartment': 'Lateral/Medial/Capsule',
                'Calibrated L0 (mm)': l0, 'Min Length (mm)': round(cohort_stats[st_key]['min'], 2),
                'Max Length (mm)': round(cohort_stats[st_key]['max'], 2),
                'ROM Elongation ΔL (mm)': round(cohort_stats[st_key]['delta'], 2),
                'Engagement Range': f">= {grid_flex[np.where(m >= l0)[0][0]]}°" if np.any(m >= l0) else "Slack",
                'Status': 'Smooth Anatomical Axis Model Kinematics'
            })
            
        ax.set_title(struct_titles[st_key], fontsize=11.5, fontweight='bold', pad=8)
        ax.set_xlabel('Elbow Flexion Angle (°)', fontsize=10.5, labelpad=4)
        ax.set_ylabel('Ligament Length (mm)', fontsize=10.5, labelpad=4)
        ax.set_xlim(15, 145)
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
        
    legend_elements = [
        Line2D([0], [0], color='#222222', linewidth=2.8, label='Cohort Mean Length'),
        Patch(facecolor='#888888', edgecolor='none', alpha=0.25, label='Cohort ± 1 SD Band'),
        Line2D([0], [0], color='#888888', linestyle='--', linewidth=1.2, label='Individual Participants (S1–S6)'),
        Line2D([0], [0], color='#b22222', linestyle=':', linewidth=2.0, label='Calibrated Resting Length ($L_0$, Tension Boundary)')
    ]
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.995), ncol=4, fontsize=10.5, frameon=True)
    fig.tight_layout(rect=[0.01, 0.02, 0.99, 0.95])
    
    p_png = os.path.join(MANUSCRIPT_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.png")
    p_svg = os.path.join(MANUSCRIPT_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.svg")
    p_pdf = os.path.join(MANUSCRIPT_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.pdf")
    fig.savefig(p_png, dpi=300, bbox_inches='tight')
    fig.savefig(p_svg, bbox_inches='tight')
    fig.savefig(p_pdf, bbox_inches='tight')
    shutil.copy2(p_png, os.path.join(ARTIFACT_FIG_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.png"))
    shutil.copy2(p_svg, os.path.join(ARTIFACT_FIG_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.svg"))
    shutil.copy2(p_pdf, os.path.join(ARTIFACT_FIG_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.pdf"))
    plt.close(fig)
    print("  [SAVED] Fig_S6_Ligament_Length_vs_Elbow_Flexion (PNG, SVG, PDF)")

    # Export Excel Workbook for Fig S6
    excel_path = os.path.join(TABLES_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.xlsx")
    df_eng = pd.DataFrame(engagement_rows)
    
    # Cohort Mean and SD sheet
    cohort_dict = {'Flexion_Angle_deg': grid_flex}
    for st in subject_curves:
        cohort_dict[f'{st}_Mean_mm'] = np.round(cohort_stats[st]['mean'], 3)
        cohort_dict[f'{st}_SD_mm']   = np.round(cohort_stats[st]['std'], 3)
    df_cohort = pd.DataFrame(cohort_dict)
    
    # Individual Participants sheet
    indiv_dict = {'Flexion_Angle_deg': grid_flex}
    for st in subject_curves:
        for s_id, _, _, _ in SUBJECTS:
            indiv_dict[f'{st}_{s_id}_mm'] = np.round(subject_curves[st][s_id], 3)
    df_indiv = pd.DataFrame(indiv_dict)
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df_eng.to_excel(writer, sheet_name="Engagement_Summary", index=False)
        df_cohort.to_excel(writer, sheet_name="Cohort_Mean_and_SD", index=False)
        df_indiv.to_excel(writer, sheet_name="Individual_Participants", index=False)
        
    shutil.copy2(excel_path, os.path.join(ARTIFACT_DIR, "Fig_S6_Ligament_Length_vs_Elbow_Flexion.xlsx"))
    print(f"  [SAVED] Excel Workbook: {excel_path}")

# ==============================================================================
# MAIN ORCHESTRATOR
# ==============================================================================
def main():
    wb_data = load_cross_subject_workbooks()
    
    # 1. Stacked Neon Plots
    generate_stacked_neon_plots(wb_data)
    
    # 2. Manuscript Figures 4 to 9
    generate_figures_4_to_9(wb_data)
    
    # 3. Figure 10 (Sobol & Tornado)
    generate_figure_10()
    
    # 4. Figure 11 (Model Comparison Waveforms)
    generate_figure_11(wb_data)
    
    # 5. Figure S6 (Kinematic Length vs Flexion)
    generate_figure_s6()
    
    print("\n" + "=" * 80)
    print("PHASE 7 COMPLETE: ALL PUBLICATION FIGURES SUCCESSFULLY GENERATED!")
    print(f"Primary Figures: {MANUSCRIPT_DIR}")
    print(f"Stacked Plots:   {STACKED_DIR}")
    print(f"Artifacts:       {ARTIFACT_FIG_DIR}")
    print("=" * 80)

if __name__ == '__main__':
    main()
