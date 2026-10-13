import glob, os, pandas as pd

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))
rows = []
for tag, cov in [
    ("_cov0.20", "20%"),
    ("_cov0.50", "50%"),
    ("_covafford", "Maximum affordable"),
    ("_asp_cov0.20", "Aspirational 20%"),
]:
    d = pd.concat(
        pd.read_csv(f) for f in glob.glob(f"results/figure4_mc{tag}_shard*.csv")
    )
    yr = [f"y{i}" for i in range(1, 31)]
    q = d.groupby(["iso3", "arm"])[yr].quantile([0.5, 0.025, 0.975])
    q.index.names = ["iso3", "arm", "q"]
    q = q.reset_index().melt(
        id_vars=["iso3", "arm", "q"], var_name="year", value_name="v"
    )
    q["year"] = q.year.str[1:].astype(int)
    q["coverage"] = cov
    rows.append(q)
out = (
    pd.concat(rows)
    .pivot_table(index=["coverage", "iso3", "arm", "year"], columns="q", values="v")
    .reset_index()
)
out.columns = ["coverage", "iso3", "arm", "year", "lo", "med", "hi"]
out.to_csv("results/itt_1012/annual_deaths_averted.csv", index=False)
print(len(out), "rows")
