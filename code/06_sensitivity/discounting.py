import glob, numpy as np, pandas as pd

d = pd.concat(
    pd.read_csv(f) for f in glob.glob("results/figure4_mc_cov0.20_shard*.csv")
)
yr = [f"y{i}" for i in range(1, 31)]
rows = []
for rate in [0, 0.03, 0.05]:
    w = 1 / (1 + rate) ** np.arange(0, 30)
    d["disc"] = d[yr].values @ w
    m = (
        d.pivot_table(index=["iso3", "draw"], columns="arm", values="disc")
        .groupby("iso3")
        .median()
    )
    for iso, x in m.iterrows():
        rows.append(
            dict(iso3=iso, rate=rate, diet=x.diet, glp1=x.glp1, ratio=x.diet / x.glp1)
        )
t = pd.DataFrame(rows)
t.to_csv("results/itt_1012/sens_discount_itt.csv", index=False)
for rate, g in t.groupby("rate"):
    print(
        f"Discounting {rate:.0%}: global ratio {g[g.iso3=='GLB'].ratio.iloc[0]:.2f}, diet leads {(g[g.iso3!='GLB'].ratio>1).sum()}/19"
    )
