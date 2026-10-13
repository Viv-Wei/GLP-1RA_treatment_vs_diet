import pandas as pd, numpy as np, glob, os

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))
PAT = {
    "20%": "results/figure4_mc_cov0.20_shard*.csv",
    "50%": "results/figure4_mc_cov0.50_shard*.csv",
    "Maximum affordable": "results/figure4_mc_covafford_shard*.csv",
}
q = lambda x: (np.median(x), np.percentile(x, 2.5), np.percentile(x, 97.5))
rows = []
for cov, pat in PAT.items():
    files = sorted(glob.glob(pat))
    if not files:
        print("no files for", cov, pat)
        continue
    mc = pd.concat([pd.read_csv(f) for f in files])
    for iso in sorted(mc.iso3.unique()):
        s = {
            a: mc[(mc.iso3 == iso) & (mc.arm == a)].set_index("draw").cum30
            for a in ["diet", "glp1"]
        }
        idx = s["diet"].index.intersection(s["glp1"].index)
        d, g = s["diet"][idx].values, s["glp1"][idx].values
        lb = np.maximum(d, g)
        add_diet = (lb - g) / g * 100
        add_glp1 = (lb - d) / d * 100
        r = dict(iso3=iso, coverage=cov, n_draws=len(idx))
        for k, v in [
            ("combo_lower_bound", lb),
            ("add_diet_to_treatment_pct", add_diet),
            ("add_treatment_to_diet_pct", add_glp1),
        ]:
            m, lo, hi = q(v)
            r.update({f"{k}_med": m, f"{k}_lo": lo, f"{k}_hi": hi})
        r["pct_draws_diet_gt_treatment"] = (d > g).mean() * 100
        rows.append(r)
out = pd.DataFrame(rows)
out.to_csv("tables/tableS19_full_redundancy_bound.csv", index=False)
cols = [
    "iso3",
    "coverage",
    "add_diet_to_treatment_pct_med",
    "add_diet_to_treatment_pct_lo",
    "add_diet_to_treatment_pct_hi",
    "pct_draws_diet_gt_treatment",
]
print(out[out.iso3 == "GLB"][cols].round(1).to_string(index=False))
for cov in PAT:
    x = out[(out.coverage == cov) & (out.iso3 != "GLB")]
    print(
        f"{cov}: countries where adding diet adds >0 at median under full redundancy: {(x.add_diet_to_treatment_pct_med > 0).sum()}/19"
    )
