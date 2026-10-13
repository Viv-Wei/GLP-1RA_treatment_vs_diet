# Overview

***GLP-1 receptor agonist treatment versus dietary policy: 30-year projections in 19 countries***

The model compares deaths averted over 30 years (2023–2053) by population-wide dietary policy and by GLP-1 receptor agonist (GLP-1RA) treatment for adults with BMI ≥25, in 19 countries and globally, and estimates the treatment coverage needed for GLP-1RA treatment to avert as many deaths as dietary policy.

---

## Root scripts
-  `run_main_pipeline.sh`  Orchestrates the projection pipeline: life-table runs (deterministic and Monte Carlo), Supplementary Tables 9-10, discounting, deterministic/structural/population sensitivities, expenditure, overlap, equity and figures; writes results under `results/`.
-  `run_static_mediator.sh`  Orchestrates the mediator-level and static comparisons: per country x outcome depth runs, affordable-coverage depth, depth aggregation, and the CKD/T2D and IHD/stroke Monte Carlo shards merged into `results/figure2b_*.csv`.

## 01_burden/
-  `build_factor_channel_burden.py`  Builds diet-attributable burden by factor and mediator channel from GBD dietary-risk CSVs (`data/gbd/`, `params/country_map.csv`, `params/glp1_channels_gbd2023.csv`); writes `results/factor_burden_by_age.csv` and `results/figure1_channel_burden.csv`.
-  `build_joint_burden.py`  Combines diet-attributable deaths for six outcomes from GBD CSVs and computes crude and age-standardised rates; writes `results/joint_burden_2023.csv` and `results/supp_table2_joint.csv`.

## 02_intervention_arms/
-  `build_calibration_k.py`  Derives the calibration factor k = official / modelled PAF per country x outcome from `results/paf_check_vs_official_v2.csv`; writes `params/paf_calibration_k.csv`.
-  `build_glp1_flow_table.py`  Computes dietary PAFs and the GLP-1 reachable/mediated flow fractions; reads GBD CSVs, `results/figure1_channel_burden.csv`, `results/glp1_arm_pif.csv` and `params/overweight_prevalence.csv`; writes `results/figure2_flow_v2.csv`.
-  `build_overweight_prevalence.py`  Derives population-weighted overweight (BMI >= 25) prevalence from GBD BMI exposure, population and WPP data; writes `params/overweight_prevalence.csv`.
-  `check_paf_official.py`  Validates modelled PAFs (`results/paf_gbdexp_2023_{male,female}.csv`) against official GBD PAFs; writes `results/paf_check_vs_official.csv`.
-  `diet_policy_arm.py`  Computes base and scenario dietary PAFs and averted fractions using dose-response curves, exposure and the sodium-via-SBP pathway; writes `results/diet_policy_arm.csv`.
-  `glp1_arm.py`  Computes GLP-1 PIFs for IHD and ischaemic stroke via SBP and FPG from trial parameters and exposure; writes `results/glp1_arm_pif.csv`.

## 03_mediation/
-  `calc_depth.py`  Computes the per country x outcome offset depth (GLP-1 shift / dietary push) for SBP/FPG/LDL at 20% coverage; writes `params/figure3_depth_{iso3}_{outcome}.json`.
-  `calc_depth_050.py`  Same as `calc_depth.py` at 50% coverage; writes `params/depth_cov050/figure3_depth_{iso3}_{outcome}.json`.
-  `calc_depth_afford.py`  Computes offset depth at the maximum affordable coverage (`params/coverage_afford_cap.csv`) for countries x cardiovascular outcomes; writes `params/depth_afford/figure3_depth_{iso3}_{outcome}.json`.
-  `gen_depth_cv_multicov.py`  Segments the Figure 3 flows into offset / mediated-not-offset / no-pathway using the depth JSONs; writes `results/f3cv_flows_*.csv`, `results/f3cv_depth_*.csv`, `results/depth_cv_country.csv` and `results/depth_cv_factor.csv`.
-  `paf_deficiency_v2.py`  Core module computing age-specific deficiency PAFs (mean and distributional) from dose-response curves, exposure and the channel/convergence tables; writes `results/paf_denominator.csv`.
-  `run_mediator_figures.sh`  Runs `gen_depth_cv_multicov.py` and then renders the Sankey and factor-bar mediation figures for GLB at 20%, 50% and affordable coverage.

## 04_static_comparison/
-  `mc_static_ckd_t2d.py`  Monte Carlo (200 CKD / 100 T2D draws) of diet vs GLP-1 avertable shares for CKD and T2D; writes `results/figure2b_ckd_t2d_draws/arms_ui/ratio_ui.csv`.
-  `mc_static_ihd_is.py`  One shard of the IHD/ischaemic stroke Monte Carlo (1000 / nshards draws); writes `results/figure2b_arms_ui_19c_shard{N}.csv` and `results/f2b_draws_shard{N}.csv`.
-  `merge_static_shards.py`  Merges the 8 IHD/stroke shard CSVs into `results/figure2b_arms_ui_19c.csv` and reports between-shard variability.

## 05_projection/
-  `build_tables_9_10.py`  Aggregates the Figure 4 Monte Carlo shards into `results/itt_1012/tableS9_deaths_itt.csv` and `results/itt_1012/tableS10_ranking_threshold_itt.csv`.
-  `expenditure_per_death_averted.py`  Computes GLP-1 vs diet cost per death averted from the cumulative-runs table, population and prices; writes `tables/tableS_cost_analysis.csv`.
-  `lifetable_full.py`  Deterministic three-state life-table simulation of the diet, GLP-1 and combination arms from GBD incidence/prevalence/deaths and the flow and policy tables; writes `results/figure4_full.csv`.
-  `lifetable_mc.py`  Runs sharded Monte Carlo draws by perturbing the flow and policy inputs and calling `lifetable_full.py`; writes `results/figure4_mc{tag}_cov{cov}_shard{N}.csv`.
-  `summarise_annual_deaths.py`  Summarises annual quantiles across the Figure 4 Monte Carlo shards; writes `results/itt_1012/annual_deaths_averted.csv`.
-  `t2d_incidence_arm.py`  Projects averted new T2D cases per year for the GLP-1 and diet arms from GBD T2D incidence, the trial hazard ratio and the dietary PAF; writes `results/t2d_incidence_arms.csv`.

## 06_sensitivity/
-  `discounting.py`  Applies 0 / 3 / 5% annual discounting to the Figure 4 Monte Carlo cumulative deaths; writes `results/itt_1012/sens_discount_itt.csv`.
-  `lifetable_structural_sensitivity.py`  Reruns the life table with structural variants (same-cycle case fatality; diseased-compartment background-mortality multipliers) across countries; writes `results/itt_1012/sens_structural_countries_itt.csv` and `tables/tableS_lifetable_structural_sensitivity.csv`.
-  `overlap_full_redundancy_bound.py`  Computes full-redundancy lower bounds for the diet + GLP-1 combination from the Figure 4 Monte Carlo shards; writes `tables/tableS19_full_redundancy_bound.csv`.
-  `overlap_sensitivity.py`  Sweeps 0-100% overlap between the dietary and GLP-1 effects using the Figure 4 Monte Carlo shards; writes `tables/tableS16_overlap_sensitivity.csv`.
-  `persistence_implementation_sweep.py`  Runs a 48-scenario grid over policy maintenance, GLP-1 retention and diet ramp-up through `lifetable_full.py`; writes `results/sensitivity_sweep_19c.csv`.
-  `population_sensitivity.py`  Reruns the life table with the 2053 UN WPP population structure instead of 2023; writes `tables/tableS_population_sensitivity_v2.csv`.
-  `run_deterministic_sensitivity.py`  Runs deterministic sensitivity sub-jobs (incidence split; cardiovascular and CKD hazard ratios) through `lifetable_full.py`; writes `results/itt_1012/sens_deterministic_itt.csv`.

## 07_equity/
-  `affordability_model.py`  Merges income quintiles (`params/income_quintiles_wb.csv`) with GLP-1 prices to compute cost as a percentage of income; writes `results/affordability_matrix.csv`.
-  `build_equity_twoarm.py`  Allocates dietary and GLP-1 deaths averted across education groups using Wittgenstein shares, equity PAFs, joint burden and affordability; writes `tables/tableS_equity_twoarm.csv` and `tables/tableS_equity_mapping_sensitivity.csv`.
-  `equity_burden_by_edu.py`  Combines education PAFs, Wittgenstein population shares and the Figure 1 channel burden; writes `results/equity_burden_by_edu.csv`.
-  `equity_mc.py`  Monte Carlo of the low/high-education PAF ratio from exposure gradients and dose-response curves; writes `results/equity_ratio_ui.csv`.
-  `equity_paf_by_edu.py`  Computes education-specific PAFs by scaling exposure with the cross-education / national intake ratios; writes `results/equity_paf_by_edu.csv`.

## figures/
-  `plot_figure1_starting_points.py`  Draws Figure 1 (age-standardised rates by outcome, overweight prevalence, dietary-factor shares) from the Supplementary Tables workbook and `params/overweight_prevalence.csv`; writes `figures/figure1_starting_points.*`.
-  `plot_figure2_projections.py`  Draws the global and per-country annual/cumulative projection panels from the `results/itt_1012/` tables; writes `figures/figure2_projections*.{pdf,png}`.
-  `plot_figure3_threshold.py`  Plots the dietary-vs-GLP-1 threshold coverage curves and per-country thresholds from the `results/itt_1012/` tables and sensitivity sweep; writes `figures/figure3_threshold.*`.
-  `plot_figure4_country_drivers.py`  Relates targeted dietary share, overweight prevalence and threshold coverage and draws a factor-priority heat map; writes `figures/figure4_country_drivers.*`.
-  `plot_figure5a_robustness.py`  Compiles all sensitivity results (deterministic, structural, persistence, population, discounting) into a forest-style robustness figure; writes `figures/figure5_robustness.*`.
-  `plot_figure5b_sensitivity_heatmap.py`  Builds a country x scenario heat map of diet/GLP-1 ratios; writes `figures/suppfig_sensitivity_heatmap.*`.
-  `plot_figure6_equity.py`  Draws the equity panels (PAF ratio, affordability, education-group slope charts); writes `figures/figure6_equity*.*`.
-  `plot_suppfig2_mediation.py`  Draws the mediation-category stacked bars and the accounting-rule dot plot; writes `figures/suppfig2_mediation.{pdf,png}`.
-  `plot_suppfig3a_sankey.R`  Renders the factor-mediator-outcome alluvial/Sankey diagram per country x coverage from `results/f3cv_flows_*.csv`; writes `figures/figure3_all/sankey_cv_*.{pdf,png}`.
-  `plot_suppfig3b_factor_bars.R`  Renders the stacked dietary-factor bars by offset segment from `results/f3cv_depth_*.csv`; writes `figures/figure3_all/foodbar_cv_*.{pdf,png}`.
-  `plot_suppfig4_static.py`  Draws per-outcome static avertable-share comparisons (diet vs GLP-1 by scenario); writes `figures/suppfig4_static.{pdf,png}`.
-  `plot_suppfig8_expenditure.R`  Renders radial bar charts of diet vs GLP-1 expenditure per death averted from `tables/tableS_cost_analysis.csv`; writes `figures/suppfig8_expenditure.{pdf,png}`.

---

## Data and parameters
-  `params/`  Model parameters: country map, dose-response curves, mediation matrix, trial effects, prices, income quintiles, overweight prevalence and exposures.
-  `results/`  Final results used by the figure scripts; `results/itt_1012/` holds the projection and sensitivity outputs.

**Raw data**: Global Burden of Disease Study 2023 (IHME), Global Dietary Database, UN World Population Prospects 2024, Wittgenstein Centre Human Capital Data Explorer, World Bank income distribution data, and published trial results (SELECT, FLOW, SUSTAIN-6, SOUL, STEP 1). Download instructions are in `params/Raw_data.md`. These data are subject to their providers' terms of use and are not redistributed here.
