import pandas as pd, numpy as np, csv, glob, os

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))

rows_s16 = []
coverage_configs = [
    ("20%", "0.20", "results/figure4_mc_cov0.20_shard*.csv"),
    ("50%", "0.50", "results/figure4_mc_cov0.50_shard*.csv"),
    ("Affordable", "afford", "results/figure4_mc_covafford_shard*.csv"),
]

for cov_label, cov_id, pattern in coverage_configs:
    mc_files = sorted(glob.glob(pattern))
    if not mc_files:
        print(f"WARNING: no files for {cov_label} ({pattern})")
        continue
    mc_draws = pd.concat([pd.read_csv(f) for f in mc_files])
    print(f"Loaded {cov_label}: {len(mc_draws)} rows from {len(mc_files)} shards")

    for iso in sorted(mc_draws.iso3.unique()):
        diet_d = (
            mc_draws[(mc_draws.iso3 == iso) & (mc_draws.arm == "diet")]
            .set_index("draw")
            .cum30
        )
        glp1_d = (
            mc_draws[(mc_draws.iso3 == iso) & (mc_draws.arm == "glp1")]
            .set_index("draw")
            .cum30
        )
        combo_d = (
            mc_draws[(mc_draws.iso3 == iso) & (mc_draws.arm == "combo")]
            .set_index("draw")
            .cum30
        )
        common = diet_d.index.intersection(glp1_d.index).intersection(combo_d.index)
        if len(common) == 0:
            continue
        diet = diet_d[common]
        glp1 = glp1_d[common]
        combo = combo_d[common]

        max_overlap = np.minimum(diet * 0.72, glp1 * 0.67)
        for pct in [0, 25, 50, 75, 100]:
            combo_adj = np.maximum(
                combo - max_overlap * pct / 100, np.maximum(diet, glp1)
            )
            add_ben = (combo_adj - glp1) / glp1 * 100
            rows_s16.append(
                {
                    "Country": iso,
                    "Coverage": cov_label,
                    "Overlap (%)": pct,
                    "Combination adjusted (median)": round(combo_adj.median()),
                    "Combination adjusted (95% UI lower)": round(
                        combo_adj.quantile(0.025)
                    ),
                    "Combination adjusted (95% UI upper)": round(
                        combo_adj.quantile(0.975)
                    ),
                    "Additional benefit over GLP-1 alone, median (%)": round(
                        add_ben.median(), 1
                    ),
                    "Additional benefit over GLP-1 alone, 95% UI lower (%)": round(
                        add_ben.quantile(0.025), 1
                    ),
                    "Additional benefit over GLP-1 alone, 95% UI upper (%)": round(
                        add_ben.quantile(0.975), 1
                    ),
                }
            )

df_s16 = pd.DataFrame(rows_s16)
with open("tables/tableS16_overlap_sensitivity.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(
        [
            "Supplementary Table 16. Overlap sensitivity (0-100%) with 95% UI from MC draws."
        ]
    )
    w.writerow(
        [
            "Overlap fraction of GLP-1 unexplained (67%) assumed to duplicate dietary direct pathway (72%)."
        ]
    )
    w.writerow(
        [
            "Additional benefit = (combination adjusted - GLP-1 alone) / GLP-1 alone x 100."
        ]
    )
    w.writerow([])
    w.writerow(df_s16.columns.tolist())
    for _, r in df_s16.iterrows():
        w.writerow(r.tolist())

print(f"\nSaved: {len(df_s16)} rows")
print("\nGLB results at 100% overlap:")
glb = df_s16[(df_s16.Country == "GLB") & (df_s16["Overlap (%)"] == 100)]
for _, r in glb.iterrows():
    print(
        f"  {r['Coverage']}: +{r['Additional benefit over GLP-1 alone, median (%)']:.1f}% ({r['Additional benefit over GLP-1 alone, 95% UI lower (%)']:.1f}-{r['Additional benefit over GLP-1 alone, 95% UI upper (%)']:.1f}%)"
    )
