set -euo pipefail
cd "${CRA_ROOT:-$(pwd)}"
source config/env_main.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
TAG=${TAG:-$(date +%m%d_%H%M)}
step() { echo; echo "===== [$(date +%H:%M:%S)] $*"; }

step "A0. Back up depth files -> params/prev_depth_${TAG}/"
B=params/prev_depth_${TAG}; mkdir -p $B
cp params/figure3_depth_*.json $B/ 2>/dev/null || true
cp -r params/depth_cov050 params/depth_afford $B/ 2>/dev/null || true
cp results/depth_cv_country.csv results/depth_cv_factor.csv $B/ 2>/dev/null || true
mkdir -p params/depth_cov050 params/depth_afford


step "A1. Depth at 20% and 50% coverage (all countries x 4 cardiovascular outcomes, parallel)"
ISOS=$(ls results/figure3_triple_*.csv | grep -v depth | sed -E 's#.*figure3_triple_(.*)\.csv#\1#')
: > /tmp/depth_jobs.txt
for iso in $ISOS; do for o in IHD IschStroke ICH SAH; do
  echo "python code/03_mediation/calc_depth.py --iso3 $iso --outcome $o"                     >> /tmp/depth_jobs.txt
  echo "python code/03_mediation/calc_depth_050.py --iso3 $iso --outcome $o"     >> /tmp/depth_jobs.txt
done; done
echo "jobs: $(wc -l < /tmp/depth_jobs.txt)"
cat /tmp/depth_jobs.txt | xargs -P 64 -I{} bash -c '{}' > logs/depth_${TAG}.log 2>&1
echo "20% files: $(ls params/figure3_depth_*.json | wc -l) | 50% files: $(ls params/depth_cov050/*.json | wc -l)"

step "A2. Depth at maximum affordable coverage"
python code/03_mediation/calc_depth_afford.py > logs/depth_afford_${TAG}.log 2>&1; tail -2 logs/depth_afford_${TAG}.log

step "A3. Supplementary Tables 5-6 (offset by coverage, by country and factor)"
python code/03_mediation/gen_depth_cv_multicov.py > logs/depth_summary_${TAG}.log 2>&1; head -25 logs/depth_summary_${TAG}.log
python3 -c "
import pandas as pd
o = pd.read_csv('$B/depth_cv_country.csv'); n = pd.read_csv('results/depth_cv_country.csv')
m = o.merge(n, on=['iso3','coverage'], suffixes=('_old','_new'))
print(m[m.iso3=='GLB'][['coverage','offset_pct_old','offset_pct_new','mediated_not_offset_pct_new','no_pathway_pct_new']].to_string(index=False))"

step "B0. Back up static comparison -> results/prev_static_${TAG}/"
S=results/prev_static_${TAG}; mkdir -p $S
mv results/figure2b_arms_ui_19c_shard*.csv results/f2b_draws_shard*.csv $S/ 2>/dev/null || true
cp results/figure2b_arms_ui_19c.csv results/figure2b_ckd_t2d_arms_ui.csv results/figure2b_ckd_t2d_ratio_ui.csv $S/ 2>/dev/null || true

step "B1. CKD and T2D: 300 draws on 100 workers (parallel version of mc_figure2b_ckd_t2d_fast_v3.py)"
python code/04_static_comparison/mc_static_ckd_t2d.py --workers 100 > logs/static_ckd_t2d_${TAG}.log 2>&1
tail -8 logs/static_ckd_t2d_${TAG}.log

step "B2. IHD and IS: 64 shards in parallel, then merge"
seq 0 63 | xargs -P 64 -I{} python code/04_static_comparison/mc_static_ihd_is.py --shard {} --nshards 64 > logs/static_ihd_is_${TAG}.log 2>&1
n=$(ls results/figure2b_arms_ui_19c_shard*.csv | wc -l); echo "shards: $n / 64"
[ "$n" -eq 64 ] || { echo "missing shards - stop"; exit 1; }
python code/04_static_comparison/merge_static_shards.py

step "Done. Backups: $B and $S"
