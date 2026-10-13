set -euo pipefail
export CRA_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$CRA_ROOT"
source config/env_main.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
STAGE="${1:-all}"
mkdir -p logs figures results/itt_1012
step() { echo; echo "===== [$(date +%H:%M:%S)] $*"; }

if [[ "$STAGE" == all || "$STAGE" == params ]]; then
  step "params: overweight prevalence and GLP-1RA flow table (raw GBD and WPP files required)"
  python code/02_intervention_arms/build_overweight_prevalence.py
  python code/02_intervention_arms/build_glp1_flow_table.py
fi

if [[ "$STAGE" == all || "$STAGE" == main ]]; then
  step "main: projections, Monte Carlo, sensitivity analyses, expenditure, overlap, equity"
  bash code/run_main_pipeline.sh
  python code/05_projection/summarise_annual_deaths.py
fi

if [[ "$STAGE" == all || "$STAGE" == static ]]; then
  step "static: mediator-level and static comparisons"
  bash code/run_static_mediator.sh
fi

if [[ "$STAGE" == all || "$STAGE" == mediator_figures ]]; then
  step "mediator figures: Supplementary Figure 3 panels"
  bash code/03_mediation/run_mediator_figures.sh
fi

if [[ "$STAGE" == all || "$STAGE" == figures ]]; then
  step "figures: Figures 1-6, Supplementary Figures 2 and 4-8"
  for s in plot_figure1_starting_points plot_figure2_projections plot_figure3_threshold plot_figure4_country_drivers \
           plot_figure5a_robustness plot_figure5b_sensitivity_heatmap plot_figure6_equity \
           plot_suppfig2_mediation plot_suppfig4_static; do
    python "code/figures/$s.py"
  done
  Rscript code/figures/plot_suppfig8_expenditure.R
fi
step "done"
