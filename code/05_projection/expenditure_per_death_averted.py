import pandas as pd, sys, os

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))
sys.path.insert(0, "code/05_projection")
from lifetable_full import load_pop

p = load_pop()
isos = sorted(set(k[0] for k in p.keys()))
adults = {i: sum(v for k, v in p.items() if k[0] == i) for i in isos}
ow = dict(
    pd.read_csv("params/overweight_prevalence.csv")[["iso3", "overweight_prev"]].values
)
brand = dict(pd.read_csv("params/glp1_prices.csv")[["iso3", "brand_annual_usd"]].values)
mc = pd.read_csv("tables/figure4/TableS_cum30_all_runs.csv")
mc = mc[mc.run == "0.20"]
PERSIST, DIET_CPP, YEARS = 0.50, 0.55, 30
rows = []
for iso in isos:
    g = mc[(mc.iso3 == iso) & (mc.arm == "glp1")].med.values
    d = mc[(mc.iso3 == iso) & (mc.arm == "diet")].med.values
    if not len(g) or not len(d):
        continue
    treated = adults[iso] * ow.get(iso, 0.3) * 0.20
    cd = adults[iso] * DIET_CPP * YEARS / d[0]
    for label, price in [
        ("Branded", brand.get(iso)),
        ("Generic $140/yr", 140),
        ("Generic $28/yr", 28),
    ]:
        if price is None or pd.isna(price):
            continue
        cg = treated * price * (1 + (YEARS - 1) * PERSIST) / g[0]
        rows.append(
            {
                "Country": iso,
                "Coverage": "20%",
                "Price scenario": label,
                "GLP-1 cost/death (USD)": round(cg),
                "Diet cost/death (USD)": round(cd),
                "Ratio GLP-1/Diet": round(cg / cd, 1),
            }
        )
out = pd.DataFrame(rows)
out.to_csv("tables/tableS_cost_analysis.csv", index=False)
print(out[out.Country.isin(["GLB", "AUS"])].to_string(index=False))
b = out[out["Price scenario"] == "Branded"]
print("\nBranded > $1M:", (b["GLP-1 cost/death (USD)"] > 1e6).sum(), "of", len(b))
