import pandas as pd, glob, numpy as np

ISO = {
    "China": "CHN",
    "India": "IND",
    "Brazil": "BRA",
    "Nigeria": "NGA",
    "United States of America": "USA",
    "Global": "GLB",
}
DIS = {
    "Ischemic heart disease": "IHD",
    "Ischemic stroke": "IschStroke",
    "Diabetes mellitus type 2": "T2D",
}
AGE15 = [f"{a}-{a+4} years" for a in range(25, 95, 5)] + ["95+ years"]

off = []
for dis, tag in DIS.items():
    g = glob.glob(f"data/gbd/{dis}/Dietary risks/*/*.csv")
    if not g:
        print(f"[missing] {dis}/Dietary risks")
        continue
    d = pd.read_csv(g[0], low_memory=False)
    d = d[
        (d.metric_name == "Percent")
        & (d.measure_name == "Deaths")
        & (d.year == 2023)
        & (d.age_name.isin(AGE15))
        & (d.sex_name.isin(["Male", "Female"]))
        & (d.location_name.isin(ISO))
    ]
    d["iso3"] = d.location_name.map(ISO)
    d["outcome"] = tag
    off.append(d[["iso3", "outcome", "sex_name", "age_name", "val"]])
off = pd.concat(off).rename(columns={"val": "PAF_official"})
off["sex"] = off.sex_name.str.lower()
off["age"] = off.age_name.str.replace(" years", "").str.replace(
    "95+", "95 plus", regex=False
)

mine = []
for sex in ["male", "female"]:
    m = pd.read_csv(f"results/paf_gbdexp_2023_{sex}.csv")
    g = (
        m.groupby(["iso3", "outcome", "age"])
        .PAF_mean.apply(lambda p: 1 - np.prod(1 - p))
        .reset_index()
    )
    g["sex"] = sex
    mine.append(g)
mine = pd.concat(mine).rename(columns={"PAF_mean": "PAF_mine"})
mine["age"] = mine.age.astype(str).str.replace(" to ", "-")
off["age"] = off.age.str.replace(" to ", "-")

t = mine.merge(
    off[["iso3", "outcome", "sex", "age", "PAF_official"]],
    on=["iso3", "outcome", "sex", "age"],
    how="inner",
)
if not len(t):
    print("!! no matches; check the age format")
    print("mine:", sorted(mine.age.unique())[:4])
    print("off :", sorted(off.age.unique())[:4])
    raise SystemExit

t["diff_pp"] = (t.PAF_mine - t.PAF_official) * 100
t["ratio"] = t.PAF_mine / t.PAF_official.replace(0, np.nan)

print(f"matched {len(t):,} cells\n")
print("=== Summary by outcome x country (mean PAF %, across age and sex) ===")
s = t.groupby(["outcome", "iso3"]).agg(
    mine=("PAF_mine", lambda x: x.mean() * 100),
    official=("PAF_official", lambda x: x.mean() * 100),
    diff_pp=("diff_pp", "mean"),
    ratio=("ratio", "mean"),
)
print(s.round(2).to_string())
print("\n=== 15 cells with the largest deviation ===")
print(
    t.reindex(t.diff_pp.abs().sort_values(ascending=False).index)
    .head(15)[["iso3", "outcome", "sex", "age", "PAF_mine", "PAF_official", "diff_pp"]]
    .round(4)
    .to_string(index=False)
)
t.to_csv("results/paf_check_vs_official.csv", index=False)
print("\n-> results/paf_check_vs_official.csv")
