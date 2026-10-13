set -euo pipefail
cd "${CRA_ROOT:-$(pwd)}"
source config/env_main.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
TAG=${TAG:-$(date +%m%d_%H%M)}
P=${P:-120}
OUT=results/itt_1012
mkdir -p logs $OUT
step() { echo; echo "===== [$(date +%H:%M:%S)] $*"; }

step "0. Back up previous results -> results/prev_${TAG}/"
PREV=results/prev_${TAG}; mkdir -p $PREV/mc
mv results/figure4_mc*_shard*.csv $PREV/mc/ 2>/dev/null || true
cp -r $OUT $PREV/ 2>/dev/null || true
cp tables/*.csv $PREV/ 2>/dev/null || true
cp tables/figure4/TableS_cum30_all_runs.csv $PREV/ 2>/dev/null || true

step "1. Deterministic check (global, 20% coverage)"
python code/05_projection/lifetable_full.py --years 30 --coverage 0.20 --diet-scenario feasible --out /tmp/lt_check.csv > /dev/null 2>&1
python3 -c "
import pandas as pd; t = pd.read_csv('/tmp/lt_check.csv'); yr = [c for c in t.columns if c.startswith('yr')]
g = t.groupby(['iso3','arm'])[yr].sum().sum(axis=1).unstack('arm').loc['GLB']/1e6
print(f'diet {g.diet:.2f}M  glp1 {g.glp1:.2f}M  ratio {g.diet/g.glp1:.2f}')"

step "2. Monte Carlo: 4 scenarios x 32 shards, $P in parallel"
NS=32; : > /tmp/mc_jobs.txt
for k in $(seq 0 $((NS-1))); do
  for cov in 0.20 0.50 afford; do
    echo "python code/05_projection/lifetable_mc.py --draws 500 --shard $k --nshards $NS --coverage $cov --diet-scenario feasible" >> /tmp/mc_jobs.txt
  done
  echo "python code/05_projection/lifetable_mc.py --draws 500 --shard $k --nshards $NS --coverage 0.20 --diet-scenario aspirational" >> /tmp/mc_jobs.txt
done
cat /tmp/mc_jobs.txt | xargs -P $P -I{} bash -c '{}' > logs/mc_${TAG}.log 2>&1
n=$(ls results/figure4_mc*_shard*.csv | wc -l); echo "shards: $n / 128"
grep -c "engine failed" logs/mc_${TAG}.log && { echo "ENGINE FAILURES - stop"; exit 1; } || true
[ "$n" -eq 128 ] || { echo "missing shards - stop"; exit 1; }

step "3. Supplementary Tables 9-10, threshold coverage, discounting"
python code/05_projection/build_tables_9_10.py | tail -22
python code/06_sensitivity/discounting.py

step "4. Deterministic sensitivity analyses: split, HR, CKD excluded, 48-scenario factorial, structure, population"
python code/06_sensitivity/run_deterministic_sensitivity.py           > logs/sens_${TAG}.log 2>&1 &
python code/06_sensitivity/persistence_implementation_sweep.py      > logs/sweep48_${TAG}.log 2>&1 &
python code/06_sensitivity/lifetable_structural_sensitivity.py > logs/structural_${TAG}.log 2>&1 &
python code/06_sensitivity/population_sensitivity.py        > logs/population_${TAG}.log 2>&1 &
SENS_ELIG_FILE=params/elig_scale_bmi30_gbd2023.csv \
  python code/05_projection/lifetable_full.py --years 30 --coverage 0.20 --diet-scenario feasible --out /tmp/sens_bmi30.csv > /dev/null 2>&1 &
python code/05_projection/t2d_incidence_arm.py --coverage 0.20 --diet-scenario feasible --years 30 \
  --ramp-glp1 1 --ramp-diet 8 --out results/t2d_incidence_arms_v2.csv > logs/t2d_${TAG}.log 2>&1 &
wait
python3 -c "
import pandas as pd; t = pd.read_csv('/tmp/sens_bmi30.csv'); yr = [c for c in t.columns if c.startswith('yr')]
s = t.groupby(['iso3','arm'])[yr].sum().sum(axis=1).unstack('arm'); s['ratio'] = s.diet/s.glp1
s.to_csv('$OUT/sens_bmi30_itt.csv'); print(f\"BMI>=30: global ratio {s.loc['GLB','ratio']:.2f}, diet leads {(s.drop('GLB').ratio>1).sum()}/19\")"
python code/05_projection/t2d_mc.py > logs/t2d_mc_${TAG}.log 2>&1; tail -8 logs/t2d_mc_${TAG}.log
cat logs/sens_${TAG}.log; tail -4 logs/sweep48_${TAG}.log; tail -9 logs/structural_${TAG}.log; tail -3 logs/population_${TAG}.log

step "5. Expenditure (Table 17) and overlap (Table 24)"
python3 << 'PYEOF'
import glob, pandas as pd
F = "tables/figure4/TableS_cum30_all_runs.csv"
TAG = {"0.20": "_cov0.20", "0.50": "_cov0.50", "afford": "_covafford"}
rows = []
for run, tag in TAG.items():
    d = pd.concat(pd.read_csv(f) for f in glob.glob(f"results/figure4_mc{tag}_shard*.csv"))
    for (iso, arm), g in d.groupby(["iso3", "arm"]):
        rows.append(dict(iso3=iso, arm=arm, med=g.cum30.median(), lo=g.cum30.quantile(.025), hi=g.cum30.quantile(.975), run=run))
pd.DataFrame(rows).to_csv(F, index=False); print("rebuilt", F)
PYEOF
python code/05_projection/expenditure_per_death_averted.py
python code/06_sensitivity/overlap_full_redundancy_bound.py
python code/06_sensitivity/overlap_sensitivity.py

step "6. Equity: allocation by education group (Wittgenstein education structure, rank-based income link)"
python code/07_equity/build_equity_twoarm.py

step "7. Figures"
python code/figures/plot_figure3_threshold.py
python code/figures/plot_figure5a_robustness.py
python code/figures/plot_figure5b_sensitivity_heatmap.py
step "Done. Previous run backed up in $PREV"
