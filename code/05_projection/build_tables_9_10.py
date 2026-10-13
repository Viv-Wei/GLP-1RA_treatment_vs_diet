import glob, pandas as pd, numpy as np


def load(tag):
    return pd.concat(
        pd.read_csv(f) for f in glob.glob(f"results/figure4_mc{tag}_shard*.csv")
    )


ARM = {
    "diet": "Dietary policy (feasible)",
    "glp1": "GLP-1RA treatment",
    "combo": "Combination",
}
COV = {
    "_cov0.20": "20% coverage",
    "_cov0.50": "50% coverage",
    "_covafford": "Maximum affordable coverage",
}
rows = []
for tag, cov in COV.items():
    d = load(tag)
    for (iso, arm), g in d.groupby(["iso3", "arm"]):
        rows.append(
            dict(
                iso3=iso,
                arm=ARM[arm],
                coverage=cov,
                med=g.cum30.median(),
                lo=g.cum30.quantile(0.025),
                hi=g.cum30.quantile(0.975),
            )
        )
asp = load("_asp_cov0.20")
for iso, g in asp[asp.arm == "diet"].groupby("iso3"):
    rows.append(
        dict(
            iso3=iso,
            arm="Dietary policy (aspirational)",
            coverage="20% coverage",
            med=g.cum30.median(),
            lo=g.cum30.quantile(0.025),
            hi=g.cum30.quantile(0.975),
        )
    )
t9 = pd.DataFrame(rows).round(0)
t9.to_csv("results/itt_1012/tableS9_deaths_itt.csv", index=False)

d = load("_cov0.20")
yr = [f"y{i}" for i in range(1, 31)]
out = []
for iso, g in d.groupby("iso3"):
    di = g[g.arm == "diet"].set_index("draw")
    gl = g[g.arm == "glp1"].set_index("draw")
    di, gl = (
        di.loc[di.index.intersection(gl.index)],
        gl.loc[di.index.intersection(gl.index)],
    )
    ahead = (di.cum30 > gl.cum30).mean() * 100
    cross = [
        next(
            (
                i + 1
                for i, (a, b) in enumerate(zip(di.loc[k, yr], gl.loc[k, yr]))
                if a > b
            ),
            np.nan,
        )
        for k in di.index
    ]
    thr = 20 * di.cum30 / gl.cum30
    out.append(
        dict(
            iso3=iso,
            diet_med=di.cum30.median(),
            glp1_med=gl.cum30.median(),
            higher=(
                "Dietary policy"
                if di.cum30.median() > gl.cum30.median()
                else "GLP-1RA treatment"
            ),
            diet_ahead_pct=round(ahead, 1),
            crossover_med=np.nanmedian(cross) if not np.isnan(cross).all() else np.nan,
            threshold_med=thr.median(),
            threshold_lo=thr.quantile(0.025),
            threshold_hi=thr.quantile(0.975),
        )
    )
t10 = pd.DataFrame(out).round(1)
t10.to_csv("results/itt_1012/tableS10_ranking_threshold_itt.csv", index=False)
pd.set_option("display.width", 200)
print(t10.sort_values("threshold_med").to_string(index=False))
