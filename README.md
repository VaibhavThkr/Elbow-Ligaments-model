# Subject-Specific Elbow Ligament Model

## Anatomical Axis, 3-Fiber RCL Fan & DeVita NLSR Dynamics

This repository contains the complete OpenSim biomechanical models, subject kinematic simulation data, analytical pipelines, and reproduction scripts for the **Anatomical Axis Elbow Ligament Model**.

The model incorporates an anatomically verified flexion-extension axis and models the radial collateral ligament complex as a validated **3-fiber fan** (anterior band `aRCL`, central band `cRCL`, and posterior band `pRCL`) alongside a recalibrated lateral ulnar collateral ligament (`LUCL`), modeled via the DeVita non-linear strain-rate dependent (NLSR) constitutive law.

---

## 📁 Repository Structure

```text
Elbow-Ligaments-model/
├── README.md                           # Documentation & reproduction walkthrough
├── requirements.txt                    # Python environment dependencies
├── .gitignore                          # Git tracking exclusions
│
├── models/                             # OpenSim model files & geometries
│   ├── calibrated_devita_nlsr/         # 6 Subject-specific OpenSim models (S1–S6)
│   │   ├── Elbow_model-S1.osim         # Subject S1 calibrated model
│   │   ├── Elbow_model-S2.osim ... S6.osim
│   └── Geometry/                       # 108 Bone and cartilage mesh surface files (.stl / .vtp)
│
├── simulations/                        # OpenSim kinematic simulation trials & participant data
│   ├── flexion_sweep_0_to_150.mot      # Continuous 0° to 150° elbow sweep motion profile
│   ├── 1_S1 Bicep/                     # Subject S1 kinematic analysis outputs & Excel databases
│   ├── 2_S2 Bicep/ ... 6_S6 Bicep/     # Subjects S2 through S6
│
└── scripts/                            # Complete reproducible analytical pipeline
    ├── 01_calibration/                 # Subject-specific geometry scaling & parameter calibration
    │   ├── 01_execute_phase1_model_update.py            # Segmental scaling & resting length calibration
    │   ├── 02_solve_viscosity_eta.py                    # Calibrates dynamic viscosity parameter (eta*)
    │   └── 03_compute_subject_kinematic_excursions.py   # Computes subject-specific length excursions
    │
    ├── 02_simulation/                  # Batch dynamic simulations & length processing
    │   ├── 01_run_batch_simulations.py                  # Batch OpenSim MuscleAnalysis execution
    │   ├── 02_generate_excels.py                        # Generates subject & cross-subject workbooks
    │   ├── 03_export_all_repetition_averaged_lengths.py # Exports cycle-normalized trajectories
    │   └── 04_run_simulation_rerun.py                   # Master rerun pipeline with calibrated inputs
    │
    ├── 03_benchmarks/                  # Constitutive model comparison
    │   └── 01_run_blankevoort_comparison.py             # Elastic, Damped & NLSR 3-model comparator
    │
    ├── 04_statistics/                  # Inferential statistics & Table 4 compilation
    │   ├── 01_run_stats.py                              # Two-way RM ANOVA and Holm post-hoc routines
    │   ├── 02_generate_table4.py                        # Compiles Table 4 participant and cohort sheets
    │   └── 03_update_table4_cohort_sheets.py            # Updates summary and numeric worksheets
    │
    ├── 05_sensitivity_analysis/        # Parameter sensitivity analysis
    │   └── 01_run_gsa.py                                # Saltelli Sobol Monte Carlo variance decomposition
    │
    ├── 06_plotting/                    # Waveforms & publication figures
    │   ├── 01_generate_figures.py                       # Renders Figs 4–11 & Fig S6 (.png, .pdf, .svg)
    │   ├── 02_generate_stacked_plots_exact_style.py     # Generates 34 stacked neon plots (569 × 1144 pt)
    │   └── 03_generate_paired_comparison_plots.py       # Generates paired group multi-bundle overlays
    │
    ├── 07_audit_and_verification/      # Quality assurance & concordance audits
    │   ├── 01_audit_all_data.py                         # 10-area comprehensive data audit
    │   ├── 02_verify_concordance.py                     # Cross-checks model properties and peak forces
    │   ├── 03_cross_check_plots_vs_tables.py            # Verifies plot Y-limits vs data sheets
    │   ├── 04_verify_complete_recalibration_concordance.py # Full concordance validation
    │   └── 05_execute_all_final_deliverables.py         # End-to-end master deliverable execution
    │
    └── utilities/                      # Diagnostic and inspection scripts
        ├── check_arm_markers.py
        ├── check_current_lengths.py
        ├── check_current_models.py
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
git clone https://github.com/VaibhavThkr/Elbow-Ligaments-model.git
cd Elbow-Ligaments-model
pip install -r requirements.txt
```

### Reproducing the Pipeline

1. **Verify Complete Concordance:**
   ```bash
   python scripts/07_audit_and_verification/01_audit_all_data.py
   ```
2. **Recompute Statistics & Table 4:**
   ```bash
   python scripts/04_statistics/01_run_stats.py
   python scripts/04_statistics/02_generate_table4.py
   ```
3. **Generate Publication Figures (Figs 4–11):**
   ```bash
   python scripts/06_plotting/01_generate_figures.py
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
