# Path B: Anatomical Axis Elbow Model

## 3-Fiber RCL Fan & DeVita NLSR Ligament Dynamics Release

This repository contains the complete OpenSim biomechanical models, subject kinematic simulation data, analytical pipelines, and publication deliverables for the **Path B Anatomical Axis Model**.

Path B incorporates an anatomically verified flexion-extension axis and replaces the lumped radial collateral ligament with a validated **3-fiber fan** (anterior band `aRCL`, central band `cRCL`, and posterior band `pRCL`) alongside a recalibrated lateral ulnar collateral ligament (`LUCL`), modeled via the DeVita non-linear strain-rate dependent (NLSR) constitutive law.

---

## 📁 Repository Structure

```text
03_Path_B_Anatomical_Axis_Model/
├── README.md                           # This documentation
├── requirements.txt                    # Python environment dependencies
├── .gitignore                          # Git tracking exclusions
│
├── models/                             # OpenSim model files & geometries
│   ├── calibrated_devita_nlsr/         # 6 Subject-specific OpenSim models (S1–S6)
│   │   ├── Elbow_model-S1.osim         # Subject S1 calibrated model
│   │   ├── Elbow_model-S2.osim         # Subject S2 calibrated model
│   │   ├── Elbow_model-S3.osim         # Subject S3 calibrated model
│   │   ├── Elbow_model-S4.osim         # Subject S4 calibrated model
│   │   ├── Elbow_model-S5.osim         # Subject S5 calibrated model
│   │   └── Elbow_model-S6.osim         # Subject S6 calibrated model
│   └── Geometry/                       # 109 Bone and cartilage mesh surface files (.stl / .vtp)
│
├── simulations/                        # OpenSim kinematic simulation trials
│   ├── flexion_sweep_0_to_150.mot      # Continuous 0° to 150° elbow sweep motion profile
│   ├── 1_S1 Bicep/                     # Subject S1 kinematic analysis outputs
│   ├── 2_S2 Bicep/                     # Subject S2 kinematic analysis outputs
│   ├── 3_S3 Bicep/                     # Subject S3 kinematic analysis outputs
│   ├── 4_S4 Bicep/                     # Subject S4 kinematic analysis outputs
│   ├── 5_S5 Bicep/                     # Subject S5 kinematic analysis outputs
│   └── 6_S6 Bicep/                     # Subject S6 kinematic analysis outputs
│
├── results/                            # Analytical outputs, workbooks & publication figures
│   ├── Cross_Subject_Average/          # Cohort grand average time-series workbooks (101-point cycle)
│   │   ├── Bicep_CrossSubject_Forces.xlsx      # Mean ± SD forces across all bundles and speeds
│   │   ├── Bicep_CrossSubject_Lengths.xlsx     # Kinematic length trajectories
│   │   ├── Bicep_CrossSubject_Strains.xlsx     # Dynamic engineering strain waveforms
│   │   └── Bicep_CrossSubject_StrainRates.xlsx # Strain rate waveforms
│   │
│   ├── figures/                        # Publication vector graphics & raster renders
│   │   ├── manuscript_figures/         # Figs 4–11 & Fig S6 with pRCL integrated (.png, .pdf, .svg, .TIF)
│   │   ├── stacked_neon_plots/         # 34 Calibrated stacked neon plots (569 × 1144 pt)
│   │   │   ├── screen_569x1144/        # 1:1 Screen preview PNGs
│   │   │   ├── paired_groups/          # Multi-bundle paired group comparison overlays
│   │   │   └── lateral_focused_matched/# Dedicated focused lateral collateral pair set
│   │   └── verification/               # Phase 1 & 2 model verification plots
│   │
│   ├── tables/                         # Publication summary tables & inferential statistics
│   │   ├── Table4_Participant_and_Cohort_Peaks.xlsx         # Participant peak forces & cohort means
│   │   ├── Model_Comparison_ANOVA_MethodB_Corrected.xlsx    # Two-way RM-ANOVA (Model × Tempo)
│   │   ├── Model_x_Tempo_Values_pRCL_LUCL.xlsx              # Cell means, post-hocs, & F-statistics
│   │   ├── Cohort_Statistical_Analysis_Report.xlsx          # Descriptive statistics across trials
│   │   ├── NLSR_Repeated_Measures_ANOVA_Report.xlsx         # Greenhouse-Geisser corrected ANOVA
│   │   ├── Global_Sensitivity_Analysis_Sobol.xlsx           # Saltelli Sobol indices (S1, ST)
│   │   ├── PathB_All_Participant_Peaks_Rerun.xlsx           # Complete 108-trial peak extraction
│   │   ├── PathB_Repetition_Averaged_Lengths_and_Durations.xlsx # Real-time un-normalized kinematics
│   │   ├── Phase2_Calibrated_Ligament_Parameters.xlsx       # Subject L0, eta*, and stiffness values
│   │   ├── Phase2_RCL_3Fiber_Fan_Summary.xlsx               # aRCL, cRCL, pRCL geometry & excursion
│   │   ├── Supplementary_Table_S4_Sub_Bundles.xlsx          # Detailed sub-bundle mechanical breakdown
│   │   └── Table6_MethodA_Model_Comparison_Benchmark.xlsx   # Benchmark comparator against Method A
│   │
│   └── calibration_checkpoints/        # Historical calibration sweeps and verification audits
│       ├── phase1_verification/        # Anatomical axis alignment checks
│       └── phase2_calibration/         # Excursion parameter sweeps across S1–S6
│
└── scripts/                            # Modular computational execution pipeline
    ├── 01_calibration/                 # Model geometric update & parameter tuning
    │   ├── 01_execute_phase1_model_update.py            # Updates OpenSim models with 3-fiber RCL & LUCL
    │   ├── 02_solve_viscosity_eta.py                    # Calibrates dynamic viscosity eta*
    │   └── 03_compute_subject_kinematic_excursions.py   # Computes subject excursions and L0 bounds
    │
    ├── 02_simulation/                  # Simulation batch execution & data extraction
    │   ├── 01_run_batch_simulations_pathb.py            # Batch OpenSim MuscleAnalysis execution
    │   ├── 02_generate_excels_pathb.py                  # Generates subject & cross-subject workbooks
    │   ├── 03_export_all_repetition_averaged_lengths.py # Real-time length & phase duration calculations
    │   └── 04_run_pathb_final_rerun.py                  # Master rerun pipeline with calibrated inputs
    │
    ├── 03_benchmarks/                  # Constitutive model benchmarks
    │   └── 01_run_blankevoort_comparison_pathb.py       # Elastic, Damped & NLSR 3-model comparator
    │
    ├── 04_statistics/                  # Inferential & cohort peak statistics
    │   ├── 01_run_stats_pathb.py                        # Two-way RM ANOVA and Holm post-hoc routines
    │   ├── 02_generate_table4_pathb.py                  # Compiles Table 4 participant and cohort sheets
    │   └── 03_update_table4_cohort_sheets.py            # Cohort mean updates & concordance verification
    │
    ├── 05_sensitivity_analysis/        # Global parameter sensitivity analysis
    │   └── 01_run_gsa_pathb.py                          # Saltelli Sobol Monte Carlo variance decomposition
    │
    ├── 06_plotting/                    # Manuscript and waveform plotting routines
    │   ├── 01_generate_figures_pathb.py                 # Renders Figs 4–11 & Fig S6 (.png, .pdf, .svg)
    │   ├── 02_generate_stacked_plots_exact_style.py     # Generates 34 stacked neon plots (569 × 1144 pt)
    │   └── 03_generate_paired_comparison_plots.py       # Generates paired group multi-bundle overlays
    │
    ├── 07_audit_and_verification/      # Quality assurance & concordance audits
    │   ├── 01_audit_all_data_pathb.py                   # 10-area comprehensive data audit
    │   ├── 02_verify_concordance_pathb.py               # Cross-checks model properties and peak forces
    │   ├── 03_cross_check_plots_vs_tables.py            # Verifies plot Y-limits vs data sheets
    │   ├── 04_verify_complete_recalibration_concordance.py # Full concordance with user input checks
    │   └── 05_execute_all_final_deliverables.py         # End-to-end master deliverable execution
    │
    └── utilities/                      # Diagnostic and inspection scripts
        ├── check_arm_markers.py
        ├── check_current_lengths.py
        ├── check_generic_scaling_ratio.py
        ├── check_model_scaling.py
        ├── check_ucl_scaling.py
        ├── compare_scaling.py
        └── inspect_all_plot_ranges.py
```

---

## 🔬 Biomechanical Specifications

### 1. Lateral Collateral Ligament Complex
* **aRCL (Anterior Band):** Insertion on anterior radial head; slack in baseline flexion ($L_0$ selected to prevent unphysiological compressive tension).
* **cRCL (Central Band):** Isometric guiding bundle; near-zero length excursion across flexion ($0.5\%$ strain, peak elastic force $0.04\,\text{N}$).
* **pRCL (Posterior Band):** Primary functional restraint of the radial collateral ligament; reaches exactly $5.0\%$ strain ($3.5\text{--}4.7\,\text{N}$) at peak elbow flexion across all participants.
* **LUCL (Lateral Ulnar Collateral Ligament):** Recalibrated anatomical origin; displays physiological slackening in extension and rises to $3.0\text{--}4.3\%$ peak strain during flexion.

### 2. Calibrated Viscosity ($\eta^*$)
Calibrated against subject-specific reference cycles:
$$\eta^*_{\text{RCL}} = 5.17\,\text{MPa}\cdot\text{s} \quad (\text{aRCL, cRCL, pRCL})$$
$$\eta^*_{\text{LUCL}} = 5.63\,\text{MPa}\cdot\text{s}$$

### 3. Constitutive Formulation
The ligament force is modeled via the DeVita NLSR formulation:
$$F_{\text{total}} = F_{\text{elastic}}(\varepsilon) + F_{\text{viscous}}(\dot{\varepsilon})$$
where:
* $F_{\text{elastic}}(\varepsilon) = K \cdot \text{PCSA} \cdot \left[ 1 - \exp\left( -(\varepsilon / \beta)^\alpha \right) \right]$ (Weibull fibril recruitment)
* $F_{\text{viscous}}(\dot{\varepsilon}) = c \cdot \dot{\varepsilon}$ (dynamic viscous damping)

---

## 🚀 Getting Started

### Prerequisites
* Python 3.9+
* OpenSim 4.0+ (optional, required only for re-running batch dynamic simulations from raw `.mot` files)

### Installation
Clone the repository and install the required Python packages:
```bash
git clone https://github.com/your-username/elbow-joint-pathb.git
cd elbow-joint-pathb
pip install -r requirements.txt
```

### Reproducing the Pipeline

1. **Verify Complete Concordance:**
   ```bash
   python scripts/07_audit_and_verification/01_audit_all_data_pathb.py
   ```
2. **Recompute Statistics & Table 4:**
   ```bash
   python scripts/04_statistics/01_run_stats_pathb.py
   python scripts/04_statistics/02_generate_table4_pathb.py
   ```
3. **Generate Publication Figures (Figs 4–11):**
   ```bash
   python scripts/06_plotting/01_generate_figures_pathb.py
   ```
4. **Generate Stacked Neon Plots ($569 \times 1144\text{ pt}$):**
   ```bash
   python scripts/06_plotting/02_generate_stacked_plots_exact_style.py
   ```

---

## 📊 Harmonized Plot Groupings

All 3-panel stacked plots follow the IEEE publication canvas ratio ($569 \times 1144\text{ pt}$, 7 outward ticks per Y-axis) with matched scales across functional groups:

1. **Collateral Complex (`pRCL`, `LUCL`, `pUCL`):**
   * Force: $[0, 12]\,\text{N}$
   * Strain: $[-0.50, 0.10]$
   * Strain Rate: $[-1.20, 1.20]\,\text{s}^{-1}$
2. **Anterior/Posterior UCL Complex (`aaUCL`, `caUCL`, `paUCL`):**
   * Force: $[0, 18]\,\text{N}$
   * Strain: $[-0.50, 0.10]$
   * Strain Rate: $[-1.20, 1.20]\,\text{s}^{-1}$
3. **Posterior Capsule Pair (`plCAP`, `pmCAP`):**
   * Force: $[0, 6]\,\text{N}$
   * Strain: $[-0.75, 0.15]$
   * Strain Rate: $[-1.50, 1.50]\,\text{s}^{-1}$

---

## 📜 Citation & License

This codebase and associated biomechanical models are released for academic and research reproducibility under the MIT License.
