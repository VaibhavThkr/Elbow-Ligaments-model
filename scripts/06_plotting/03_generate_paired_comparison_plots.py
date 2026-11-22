"""
Generate Paired Comparison Stacked Neon Plots (569 x 1144 pt)
Groups:
  1. pRCL, LUCL, pUCL
  2. aaUCL, caUCL, paUCL
  3. plCAP, pmCAP
"""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from PIL import Image

BASE_DIR = r"BASE_DIR"
PROJECT_DIR = os.path.join(BASE_DIR, "PROJECT_DIR")
DATA_DIR = os.path.join(PROJECT_DIR, "results", "Cross_Subject_Average")
OUT_DIR = os.path.join(PROJECT_DIR, "results", "figures", "stacked_neon_plots", "paired_groups")
ARTIFACT_DIR = r"C:/Users/vaibh/.gemini/antigravity/brain/1480b829-4eb3-4948-8aae-83a7b237bdf1"
ARTIFACT_OUT_DIR = os.path.join(ARTIFACT_DIR, "calibrated_stacked_plots_model", "paired_groups")

for d in [OUT_DIR, ARTIFACT_OUT_DIR, os.path.join(OUT_DIR, "screen_569x1144"), os.path.join(ARTIFACT_OUT_DIR, "screen_569x1144")]:
    os.makedirs(d, exist_ok=True)

GROUPS = {
    'Collateral_Trio_pRCL_LUCL_pUCL': {
        'ligs': ['pRCL', 'LUCL', 'pUCL'],
        'colors': {'pRCL': '#ff00ff', 'LUCL': '#00ffff', 'pUCL': '#ffcc00'},
        'limits': {
            'Forces':      {'lim': (0, 12),       'fmt': '%d'},
            'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
            'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
        }
    },
    'UCL_Trio_aaUCL_caUCL_paUCL': {
        'ligs': ['aaUCL', 'caUCL', 'paUCL'],
        'colors': {'aaUCL': '#ff00ff', 'caUCL': '#00ffff', 'paUCL': '#ffcc00'},
        'limits': {
            'Forces':      {'lim': (0, 18),       'fmt': '%d'},
            'Strains':     {'lim': (-0.50, 0.10), 'fmt': '%.2f'},
            'StrainRates': {'lim': (-1.20, 1.20), 'fmt': '%.2f'}
        }
    },
    'Posterior_Capsule_Pair_plCAP_pmCAP': {
        'ligs': ['plCAP', 'pmCAP'],
        'colors': {'plCAP': '#ff00ff', 'pmCAP': '#00ffff'},
        'limits': {
            'Forces':      {'lim': (0, 6),        'fmt': '%d'},
            'Strains':     {'lim': (-0.75, 0.15), 'fmt': '%.2f'},
            'StrainRates': {'lim': (-1.80, 1.80), 'fmt': '%.2f'}
        }
    }
}

VARIABLES = ['Forces', 'Strains', 'StrainRates']
YLABELS = {
    'Forces': "Force (N)",
    'Strains': "Strain",
    'StrainRates': "Strain Rate (1/s)"
}

CONDITIONS = ['nor', 'rest']
COND_NAMES = ['Normal', 'Resting']

def set_formatting():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
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

def main():
    set_formatting()
    # Load data
    wb_data = {}
    for var in VARIABLES:
        p = os.path.join(DATA_DIR, f"Bicep_CrossSubject_{var}.xlsx")
        xls = pd.ExcelFile(p)
        wb_data[var] = {s: pd.read_excel(xls, s) for s in xls.sheet_names}

    screen_ws = os.path.join(OUT_DIR, "screen_569x1144")
    screen_art = os.path.join(ARTIFACT_OUT_DIR, "screen_569x1144")

    for gname, gcfg in GROUPS.items():
        ligs = gcfg['ligs']
        colors = gcfg['colors']
        limits = gcfg['limits']

        for cond, cond_name in zip(CONDITIONS, COND_NAMES):
            fig, axes = plt.subplots(3, 1, figsize=(8, 16))
            fig.patch.set_alpha(0.0)

            for i, var in enumerate(VARIABLES):
                ax = axes[i]
                ax.patch.set_alpha(0.0)

                for lig in ligs:
                    df = wb_data[var][lig]
                    cycle_pct = df['Cycle %'].values
                    # Plot Fast condition for each bundle
                    mean_col = f"{cond}_fast_Mean"
                    sd_col   = f"{cond}_fast_SD"

                    if mean_col in df.columns:
                        m = df[mean_col].values
                        sd = df[sd_col].values if sd_col in df.columns else np.zeros_like(m)
                        valid = ~np.isnan(m)

                        if var == 'Forces':
                            m_plot = np.maximum(m, 0.0)
                            lower = np.maximum(m - sd, 0.0)
                            upper = np.maximum(m + sd, 0.0)
                        else:
                            m_plot = m
                            lower = m - sd
                            upper = m + sd

                        ax.fill_between(cycle_pct[valid], lower[valid], upper[valid],
                                        color=colors[lig], alpha=0.18, edgecolor='none')
                        ax.plot(cycle_pct[valid], m_plot[valid], color=colors[lig],
                                linewidth=3.2, label=lig)

                # Set group limits
                cfg = limits[var]
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

                # Legend on top panel
                if i == 0:
                    leg = ax.legend(loc='upper right', frameon=False, fontsize=18)
                    for text in leg.get_texts():
                        text.set_fontweight('bold')

            fig.align_ylabels(axes)
            fig.tight_layout(h_pad=3.0)

            filename_base = f"{gname}_{cond_name}_Paired_Comparison"
            p_png = os.path.join(OUT_DIR, f"{filename_base}.png")
            p_svg = os.path.join(OUT_DIR, f"{filename_base}.svg")
            p_pdf = os.path.join(OUT_DIR, f"{filename_base}.pdf")

            fig.savefig(p_png, dpi=300, transparent=True, bbox_inches='tight')
            fig.savefig(p_svg, transparent=True, bbox_inches='tight')
            fig.savefig(p_pdf, transparent=True, bbox_inches='tight')

            shutil.copy2(p_png, os.path.join(ARTIFACT_OUT_DIR, f"{filename_base}.png"))
            shutil.copy2(p_svg, os.path.join(ARTIFACT_OUT_DIR, f"{filename_base}.svg"))
            shutil.copy2(p_pdf, os.path.join(ARTIFACT_OUT_DIR, f"{filename_base}.pdf"))

            im_300 = Image.open(p_png)
            im_screen = im_300.resize((569, 1144), Image.Resampling.LANCZOS)
            im_screen.save(os.path.join(screen_ws, f"{filename_base}.png"))
            im_screen.save(os.path.join(screen_art, f"{filename_base}.png"))

            print(f"Generated Paired Group: {filename_base}.png/.svg/.pdf")
            plt.close(fig)

    print("\nSUCCESS: All paired comparison figures generated!")

if __name__ == '__main__':
    main()
