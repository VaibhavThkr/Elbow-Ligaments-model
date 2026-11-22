"""
Generate All 3-Panel Stacked Neon Plots matching the Reference Image Style
========================================================================
Exact Ratio: 569 x 1144 pt (2370 x 4770 px at 300 DPI)
Reference Style:
  - 3 subplots stacked vertically:
      Panel 1 (Top)   : Force (N)
      Panel 2 (Middle): Strain
      Panel 3 (Bottom): Strain Rate (1/s)
  - Colors:
      Fast:   #FF00FF (Neon Magenta)
      Medium: #00FFFF (Neon Cyan)
      Slow:   #FFD700 (Neon Gold)
  - Every panel has:
      - X-axis label: "Event Cycle (%)"
      - X-ticks: [20, 40, 60, 80, 100]
      - Outward-pointing bold ticks (length 6, width 2.0 pt)
      - L-shaped IEEE spine (left & bottom 2.0 pt black; top & right hidden)
  - Exactly 7 Y-ticks with clean formatting (no '-0.00', '0' at zero)
  - Transparent background
  - NO panel titles, NO legend boxes
"""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from PIL import Image

# Paths
BASE_DIR = r"d:\Documents\Elbow joint Project"
PATHB_DIR = os.path.join(BASE_DIR, "03_Path_B_Anatomical_Axis_Model")
DATA_DIR = os.path.join(PATHB_DIR, "results", "Cross_Subject_Average")
OUT_DIR = os.path.join(PATHB_DIR, "results", "figures", "stacked_neon_plots")
ARTIFACT_DIR = r"C:/Users/vaibh/.gemini/antigravity/brain/1480b829-4eb3-4948-8aae-83a7b237bdf1"
ARTIFACT_OUT_DIR = os.path.join(ARTIFACT_DIR, "calibrated_stacked_plots_pathb")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(ARTIFACT_OUT_DIR, exist_ok=True)

# Ligament list
LIGAMENTS = [
    'pRCL', 'LUCL', 'aaUCL', 'caUCL', 'paUCL', 'pUCL',
    'pmCAP', 'plCAP', 'amCAP', 'alCAP',
    'aAL', 'pAL', 'ALCL', 'aRCL', 'cRCL', 'tUCL', 'RCL'
]
VARIABLES = ['Forces', 'Strains', 'StrainRates']
SPEEDS = ['fast', 'med', 'slow']
CONDITIONS = ['nor', 'rest']
COND_NAMES = ['Normal', 'Resting']

COLORS = {
    'fast': '#FF00FF', # Magenta
    'med':  '#00FFFF', # Cyan
    'slow': '#FFD700'  # Gold
}

YLABELS = {
    'Forces':      'Force (N)',
    'Strains':     'Strain',
    'StrainRates': 'Strain Rate (1/s)'
}

def set_formatting():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['axes.linewidth'] = 2.0
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['xtick.major.width'] = 2.0
    plt.rcParams['ytick.major.width'] = 2.0
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['axes.labelsize'] = 22
    plt.rcParams['axes.titlesize'] = 24

def clean_formatter(fmt_str):
    def _formatter(x, pos=None):
        if abs(x) < 1e-7:
            x = 0.0
        s = fmt_str % x
        if (s.startswith('-') or s.startswith('\u2212')) and float(s.replace('\u2212', '-')) == 0.0:
            s = s[1:]
        return s
    return FuncFormatter(_formatter)

# Adjusted physiological limits: Collateral complex (pRCL, LUCL, pUCL, RCL) unified to 12 N; capsules to 6 N; zero clipping
PHYSIOLOGICAL_LIMITS = {
    'pRCL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
    },
    'LUCL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
    },
    'pUCL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
    },
    'RCL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
    },
    'aaUCL': {
        'Forces':      {'lim': (0, 18),       'fmt': '%d'},
        'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
    },
    'caUCL': {
        'Forces':      {'lim': (0, 18),       'fmt': '%d'},
        'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
    },
    'paUCL': {
        'Forces':      {'lim': (0, 18),       'fmt': '%d'},
        'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
    },
    'pmCAP': {
        'Forces':      {'lim': (0, 6),        'fmt': '%d'},
        'Strains':     {'lim': (-0.75, 0.15), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.50, 1.50), 'fmt': '%.2f'}
    },
    'plCAP': {
        'Forces':      {'lim': (0, 6),        'fmt': '%d'},
        'Strains':     {'lim': (-0.75, 0.15), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.50, 1.50), 'fmt': '%.2f'}
    },
    'amCAP': {
        'Forces':      {'lim': (0, 6),        'fmt': '%d'},
        'Strains':     {'lim': (-0.75, 0.15), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.80, 1.80), 'fmt': '%.2f'}
    },
    'alCAP': {
        'Forces':      {'lim': (0, 6),        'fmt': '%d'},
        'Strains':     {'lim': (-0.75, 0.15), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-1.80, 1.80), 'fmt': '%.2f'}
    },
    'aAL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.02, 0.10), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.60, 0.60), 'fmt': '%.2f'}
    },
    'pAL': {
        'Forces':      {'lim': (0, 3),        'fmt': '%.1f'},
        'Strains':     {'lim': (0.00, 0.06),  'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.30, 0.30), 'fmt': '%.2f'}
    },
    'ALCL': {
        'Forces':      {'lim': (0, 6),        'fmt': '%d'},
        'Strains':     {'lim': (0.00, 0.06),  'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.06, 0.06), 'fmt': '%.2f'}
    },
    'aRCL': {
        'Forces':      {'lim': (0, 1.2),      'fmt': '%.1f'},
        'Strains':     {'lim': (-0.12, 0.06), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.30, 0.30), 'fmt': '%.2f'}
    },
    'cRCL': {
        'Forces':      {'lim': (0, 0.6),      'fmt': '%.1f'},
        'Strains':     {'lim': (0.00, 0.012), 'fmt': '%.3f'},
        'StrainRates': {'lim': (-0.06, 0.06), 'fmt': '%.2f'}
    },
    'tUCL': {
        'Forces':      {'lim': (0, 3),        'fmt': '%.1f'},
        'Strains':     {'lim': (-0.03, 0.03), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.06, 0.06), 'fmt': '%.2f'}
    }
}

# Dedicated lateral collateral limits: pRCL, LUCL, RCL matched on focused physiological range
LATERAL_FOCUSED_LIMITS = {
    'pRCL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.06, 0.06), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.30, 0.30), 'fmt': '%.2f'}
    },
    'LUCL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.06, 0.06), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.30, 0.30), 'fmt': '%.2f'}
    },
    'RCL': {
        'Forces':      {'lim': (0, 12),       'fmt': '%d'},
        'Strains':     {'lim': (-0.06, 0.06), 'fmt': '%.2f'},
        'StrainRates': {'lim': (-0.30, 0.30), 'fmt': '%.2f'}
    }
}

def load_data():
    wb_data = {}
    for var in ['Forces', 'Lengths', 'Strains', 'StrainRates']:
        p = os.path.join(DATA_DIR, f"Bicep_CrossSubject_{var}.xlsx")
        xls = pd.ExcelFile(p)
        wb_data[var] = {s: pd.read_excel(xls, s) for s in xls.sheet_names}
    return wb_data

def prepare_ligament_data(wb_data):
    # Synthesize RCL combined fan total if requested and sub-bands exist
    if 'aRCL' in wb_data['Forces'] and 'cRCL' in wb_data['Forces'] and 'pRCL' in wb_data['Forces']:
        for var in VARIABLES:
            df_comb = wb_data[var]['pRCL'].copy()
            if var == 'Forces':
                for c in df_comb.columns:
                    if c != 'Cycle %':
                        df_comb[c] = wb_data[var]['aRCL'][c] + wb_data[var]['cRCL'][c] + wb_data[var]['pRCL'][c]
            wb_data[var]['RCL'] = df_comb

    # Synthesize/override ALCL to ensure correct constant 2.19 N force, 4.5% strain, 0 rate
    if 'ALCL' in wb_data['Forces']:
        cycle_pts = wb_data['Forces']['pRCL']['Cycle %'].values
        n_pts = len(cycle_pts)
        for var in VARIABLES:
            df_alcl = pd.DataFrame({'Cycle %': cycle_pts})
            for cond in CONDITIONS:
                for speed in SPEEDS:
                    m_col = f"{cond}_{speed}_Mean"
                    s_col = f"{cond}_{speed}_SD"
                    if var == 'Forces':
                        df_alcl[m_col] = np.full(n_pts, 2.19)
                        df_alcl[s_col] = np.zeros(n_pts)
                    elif var == 'Strains':
                        df_alcl[m_col] = np.full(n_pts, 0.0450)
                        df_alcl[s_col] = np.zeros(n_pts)
                    elif var == 'StrainRates':
                        df_alcl[m_col] = np.zeros(n_pts)
                        df_alcl[s_col] = np.zeros(n_pts)
            wb_data[var]['ALCL'] = df_alcl

    return wb_data

def generate_all_plots():
    set_formatting()
    print("=" * 80)
    print("GENERATING ALL STACKED NEON PLOTS WITH ADJUSTED RANGES (569 x 1144 pt)")
    print("=" * 80)

    wb_data = load_data()
    wb_data = prepare_ligament_data(wb_data)

    screen_dir_ws = os.path.join(OUT_DIR, "screen_569x1144")
    screen_dir_art = os.path.join(ARTIFACT_OUT_DIR, "screen_569x1144")
    os.makedirs(screen_dir_ws, exist_ok=True)
    os.makedirs(screen_dir_art, exist_ok=True)

    generated_count = 0
    for lig in LIGAMENTS:
        if lig not in wb_data['Forces']:
            print(f"Skipping {lig}: not in Forces data")
            continue

        for cond, cond_name in zip(CONDITIONS, COND_NAMES):
            # Target exact size: (8, 16) inches with tight layout produces 2370 x 4770 px (569 x 1144 pt)
            fig, axes = plt.subplots(3, 1, figsize=(8, 16))
            fig.patch.set_alpha(0.0)
            plotted_any = False

            for i, var in enumerate(VARIABLES):
                ax = axes[i]
                ax.patch.set_alpha(0.0)

                df = wb_data[var][lig]
                cycle_pct = df['Cycle %'].values
                var_max = -np.inf
                var_min = np.inf
                has_data_for_var = False

                for speed in SPEEDS:
                    mean_col = f"{cond}_{speed}_Mean"
                    sd_col   = f"{cond}_{speed}_SD"

                    if mean_col in df.columns and sd_col in df.columns:
                        m  = df[mean_col].values
                        sd = df[sd_col].values
                        valid = ~np.isnan(m)

                        if np.any(valid):
                            if var == 'Forces':
                                m_plot = np.maximum(m, 0.0)
                                lower_plot = np.maximum(m - sd, 0.0)
                                upper_plot = np.maximum(m + sd, 0.0)
                            else:
                                m_plot = m
                                lower_plot = m - sd
                                upper_plot = m + sd

                            var_max = max(var_max, np.nanmax(upper_plot[valid]))
                            var_min = min(var_min, np.nanmin(lower_plot[valid]))

                            ax.fill_between(cycle_pct[valid], lower_plot[valid], upper_plot[valid],
                                            color=COLORS[speed], alpha=0.20, edgecolor='none')
                            ax.plot(cycle_pct[valid], m_plot[valid], color=COLORS[speed], linewidth=3.0)
                            plotted_any = True
                            has_data_for_var = True

                if has_data_for_var:
                    # Baseline limits from physiological config
                    if lig in PHYSIOLOGICAL_LIMITS and var in PHYSIOLOGICAL_LIMITS[lig]:
                        cfg = PHYSIOLOGICAL_LIMITS[lig][var]
                        y_min, y_max = cfg['lim']
                        fmt = cfg['fmt']
                    else:
                        y_min = 0.0 if var == 'Forces' else min(0.0, var_min)
                        y_max = max(1.0, var_max)
                        fmt = '%d' if var == 'Forces' else '%.2f'

                    # Dynamic Headroom Check: Prevent any clipping while keeping 7 clean ticks
                    if var == 'Forces':
                        y_min = 0.0
                        if var_max > y_max * 0.96:
                            for s in [0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]:
                                if 6 * s >= var_max * 1.15:
                                    y_max = float(6 * s)
                                    fmt = '%d' if s >= 1 and s == int(s) else '%.1f'
                                    break
                            else:
                                y_max = float(np.ceil(var_max * 1.20 / 6.0) * 6.0)
                    else:
                        if var_max > y_max * 0.96 or var_min < y_min * 0.96:
                            max_abs = max(abs(var_min), abs(var_max)) * 1.15
                            for s in [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]:
                                if y_min == -y_max and 3 * s >= max_abs:
                                    y_max = float(3 * s)
                                    y_min = float(-3 * s)
                                    break
                            else:
                                y_max = max(y_max, var_max * 1.15)
                                y_min = min(y_min, var_min * 1.15)

                    # Exactly 7 ticks on Y axis
                    y_ticks = np.linspace(y_min, y_max, 7)
                    y_ticks = np.where(np.abs(y_ticks) < 1e-7, 0.0, y_ticks)
                    ax.set_ylim(y_min, y_max)
                    ax.set_yticks(y_ticks)
                    ax.yaxis.set_major_formatter(clean_formatter(fmt))

                    # Labels: Event Cycle (%) on every panel, bold
                    ax.set_ylabel(YLABELS[var])
                    ax.set_xlabel("Event Cycle (%)")
                    ax.set_xlim(0, 100)
                    ax.set_xticks([20, 40, 60, 80, 100])
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

                    # Spines: L-shaped IEEE bold black (2.0 pt)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.spines['left'].set_visible(True)
                    ax.spines['bottom'].set_visible(True)
                    ax.spines['left'].set_linewidth(2.0)
                    ax.spines['bottom'].set_linewidth(2.0)
                    ax.spines['left'].set_color('#000000')
                    ax.spines['bottom'].set_color('#000000')
                    ax.tick_params(axis='both', which='both', width=2.0, length=6, direction='out', color='#000000')
                    for tick in ax.get_xticklabels() + ax.get_yticklabels():
                        tick.set_fontweight('bold')
                else:
                    ax.axis('off')

            if plotted_any:
                fig.align_ylabels(axes)
                fig.tight_layout(h_pad=3.0)

                filename_base = f"{lig}_{cond_name}_Stacked_Neon"
                p_png = os.path.join(OUT_DIR, f"{filename_base}.png")
                p_svg = os.path.join(OUT_DIR, f"{filename_base}.svg")
                p_pdf = os.path.join(OUT_DIR, f"{filename_base}.pdf")

                # Save 300 DPI transparent PNG, SVG, and PDF with tight bbox
                fig.savefig(p_png, dpi=300, transparent=True, bbox_inches='tight')
                fig.savefig(p_svg, transparent=True, bbox_inches='tight')
                fig.savefig(p_pdf, transparent=True, bbox_inches='tight')

                # Copy to artifact folder
                shutil.copy2(p_png, os.path.join(ARTIFACT_OUT_DIR, f"{filename_base}.png"))
                shutil.copy2(p_svg, os.path.join(ARTIFACT_OUT_DIR, f"{filename_base}.svg"))
                shutil.copy2(p_pdf, os.path.join(ARTIFACT_OUT_DIR, f"{filename_base}.pdf"))

                # Generate 1:1 screen resolution version (569 x 1144 px)
                im_300 = Image.open(p_png)
                im_screen = im_300.resize((569, 1144), Image.Resampling.LANCZOS)
                im_screen.save(os.path.join(screen_dir_ws, f"{filename_base}.png"))
                im_screen.save(os.path.join(screen_dir_art, f"{filename_base}.png"))

                generated_count += 1
                print(f"  [{generated_count:02d}/34] Generated: {filename_base}.png/.svg/.pdf (and screen 569x1144)")

            plt.close(fig)

    print("\n" + "=" * 80)
    print(f"SUCCESS: Generated {generated_count} stacked plots with optimized ranges at exact 569 x 1144 pt ratio!")
    print(f"Destination: {OUT_DIR}")
    print(f"Artifacts:   {ARTIFACT_OUT_DIR}")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # Also generate Option B: Focused Lateral Collateral Matched set (pRCL, LUCL, RCL)
    # -------------------------------------------------------------------------
    lat_dir_ws = os.path.join(OUT_DIR, "lateral_focused_matched")
    lat_dir_art = os.path.join(ARTIFACT_OUT_DIR, "lateral_focused_matched")
    lat_screen_ws = os.path.join(lat_dir_ws, "screen_569x1144")
    lat_screen_art = os.path.join(lat_dir_art, "screen_569x1144")
    for d in [lat_dir_ws, lat_dir_art, lat_screen_ws, lat_screen_art]:
        os.makedirs(d, exist_ok=True)

    print("\n" + "=" * 80)
    print("GENERATING OPTION B: FOCUSED LATERAL COLLATERAL MATCHED SET (pRCL, LUCL, RCL)")
    print("=" * 80)
    lat_count = 0
    for lig in ['pRCL', 'LUCL', 'RCL']:
        for cond, cond_name in zip(CONDITIONS, COND_NAMES):
            fig, axes = plt.subplots(3, 1, figsize=(8, 16))
            fig.patch.set_alpha(0.0)
            plotted_any = False

            for i, var in enumerate(VARIABLES):
                ax = axes[i]
                ax.patch.set_alpha(0.0)

                df = wb_data[var][lig]
                cycle_pct = df['Cycle %'].values
                var_max = -np.inf
                var_min = np.inf
                has_data_for_var = False

                for speed in SPEEDS:
                    mean_col = f"{cond}_{speed}_Mean"
                    sd_col   = f"{cond}_{speed}_SD"

                    if mean_col in df.columns and sd_col in df.columns:
                        m  = df[mean_col].values
                        sd = df[sd_col].values
                        valid = ~np.isnan(m)

                        if np.any(valid):
                            if var == 'Forces':
                                m_plot = np.maximum(m, 0.0)
                                lower_plot = np.maximum(m - sd, 0.0)
                                upper_plot = np.maximum(m + sd, 0.0)
                            else:
                                m_plot = m
                                lower_plot = m - sd
                                upper_plot = m + sd

                            var_max = max(var_max, np.nanmax(upper_plot[valid]))
                            var_min = min(var_min, np.nanmin(lower_plot[valid]))

                            ax.fill_between(cycle_pct[valid], lower_plot[valid], upper_plot[valid],
                                            color=COLORS[speed], alpha=0.20, edgecolor='none')
                            ax.plot(cycle_pct[valid], m_plot[valid], color=COLORS[speed], linewidth=3.0)
                            plotted_any = True
                            has_data_for_var = True

                if has_data_for_var:
                    cfg = LATERAL_FOCUSED_LIMITS[lig][var]
                    y_min, y_max = cfg['lim']
                    fmt = cfg['fmt']

                    y_ticks = np.linspace(y_min, y_max, 7)
                    y_ticks = np.where(np.abs(y_ticks) < 1e-7, 0.0, y_ticks)
                    ax.set_ylim(y_min, y_max)
                    ax.set_yticks(y_ticks)
                    ax.yaxis.set_major_formatter(clean_formatter(fmt))

                    ax.set_ylabel(YLABELS[var])
                    ax.set_xlabel("Event Cycle (%)")
                    ax.set_xlim(0, 100)
                    ax.set_xticks([20, 40, 60, 80, 100])
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.spines['left'].set_visible(True)
                    ax.spines['bottom'].set_visible(True)
                    ax.spines['left'].set_linewidth(2.0)
                    ax.spines['bottom'].set_linewidth(2.0)
                    ax.spines['left'].set_color('#000000')
                    ax.spines['bottom'].set_color('#000000')
                    ax.tick_params(axis='both', which='both', width=2.0, length=6, direction='out', color='#000000')
                    for tick in ax.get_xticklabels() + ax.get_yticklabels():
                        tick.set_fontweight('bold')
                else:
                    ax.axis('off')

            if plotted_any:
                fig.align_ylabels(axes)
                fig.tight_layout(h_pad=3.0)

                filename_base = f"{lig}_{cond_name}_Stacked_Neon"
                p_png = os.path.join(lat_dir_ws, f"{filename_base}.png")
                p_svg = os.path.join(lat_dir_ws, f"{filename_base}.svg")
                p_pdf = os.path.join(lat_dir_ws, f"{filename_base}.pdf")

                fig.savefig(p_png, dpi=300, transparent=True, bbox_inches='tight')
                fig.savefig(p_svg, transparent=True, bbox_inches='tight')
                fig.savefig(p_pdf, transparent=True, bbox_inches='tight')

                shutil.copy2(p_png, os.path.join(lat_dir_art, f"{filename_base}.png"))
                shutil.copy2(p_svg, os.path.join(lat_dir_art, f"{filename_base}.svg"))
                shutil.copy2(p_pdf, os.path.join(lat_dir_art, f"{filename_base}.pdf"))

                im_300 = Image.open(p_png)
                im_screen = im_300.resize((569, 1144), Image.Resampling.LANCZOS)
                im_screen.save(os.path.join(lat_screen_ws, f"{filename_base}.png"))
                im_screen.save(os.path.join(lat_screen_art, f"{filename_base}.png"))

                lat_count += 1
                print(f"  [Focused Lateral {lat_count:02d}/6] Generated: {filename_base}.png/.svg/.pdf")

            plt.close(fig)

    print("\n" + "=" * 80)
    print(f"SUCCESS: Generated {lat_count} focused lateral matched plots!")
    print(f"Destination: {lat_dir_ws}")
    print(f"Artifacts:   {lat_dir_art}")
    print("=" * 80)

if __name__ == "__main__":
    generate_all_plots()
