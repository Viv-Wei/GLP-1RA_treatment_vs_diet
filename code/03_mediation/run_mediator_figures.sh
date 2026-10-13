cd "${CRA_ROOT:-$(pwd)}"
python3 code/03_mediation/gen_depth_cv_multicov.py
for cov in "20%" "50%" "affordable"; do
  Rscript code/figures/plot_suppfig3a_sankey.R GLB "$cov"
  Rscript code/figures/plot_suppfig3b_factor_bars.R GLB "$cov"
done
