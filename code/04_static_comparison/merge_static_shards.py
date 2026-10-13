import glob, pandas as pd

fs = sorted(
    f
    for f in glob.glob("results/figure2b_arms_ui_19c_shard*.csv")
    if "badtmrel" not in f and "v2run" not in f
)
assert len(fs) == 8, f"shard count {len(fs)} != 8"
d = pd.concat(
    [pd.read_csv(f).assign(shard=i) for i, f in enumerate(fs)], ignore_index=True
)
print(f"shards combined: {len(d)} rows, scenarios {sorted(d.scenario.unique())}")

g = (
    d.groupby(["scenario", "outcome", "iso3"])[["med", "lo", "hi"]]
    .mean()
    .round(2)
    .reset_index()
)
g.to_csv("results/figure2b_arms_ui_19c.csv", index=False)
print(
    f"\nmerged -> results/figure2b_arms_ui_19c.csv  ({len(g)} rows, "
    f"{g.iso3.nunique()} units x {g.outcome.nunique()} outcomes x {g.scenario.nunique()} scenarios)"
)

sd = d.groupby(["scenario", "outcome", "iso3"]).med.std()
print(
    f"Between-shard median SD: mean={sd.mean():.3f}, max={sd.max():.3f} (too large means too few draws)"
)

print("\n=== IHD: avertable share by scenario (% of diet-attributable burden) ===")
piv = g[g.outcome == "IHD"].pivot(index="iso3", columns="scenario", values="med")
cols = [
    c
    for c in ["sodium30", "feasible", "aspirational", "glp1_total", "glp1_mediated"]
    if c in piv.columns
]
print(piv[cols].to_string())
