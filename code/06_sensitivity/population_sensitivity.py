import pandas as pd, numpy as np, csv, os, sys

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))
from importlib.machinery import SourceFileLoader

lt = SourceFileLoader("lt", "code/05_projection/lifetable_full.py").load_module()

AGE15 = lt.AGE15
pop_base = lt.load_pop()
allc = lt.load_measure("queryD_All_causes_deaths.csv", "Deaths")
f2 = pd.read_csv(os.environ.get("LT_FLOW_OVERRIDE", "results/figure2_flow_v2.csv"))
f2 = f2[f2.coverage == 0.20]
dp_raw = pd.read_csv(os.environ.get("LT_DP_OVERRIDE", "results/diet_policy_arm.csv"))
dp = (
    dp_raw[dp_raw.scenario == "feasible"]
    .groupby(["outcome", "iso3"])
    .averted_frac_alldeaths.mean()
)
if dp.max() > 1:
    dp = dp / 100

_ow = pd.read_csv("params/overweight_prevalence.csv")
_ow_dict = dict(zip(_ow.iso3, _ow.overweight_prev))

OUTFILES = {
    "IHD": "queryB_ihd.csv",
    "IschStroke": "queryB_istroke.csv",
    "CKD": "queryB_ckd_3years.csv",
}

wfull = pd.read_csv(
    "data/wpp/WPP2024_PopulationByAge5GroupSex_Medium.csv",
    usecols=["Location", "ISO3_code", "Time", "AgeGrp", "PopMale", "PopFemale"],
    low_memory=False,
)
wfull = wfull[wfull.Time.isin([2023, 2053])]

AGE_MAP = {
    "25-29 years": ["25-29"],
    "30-34 years": ["30-34"],
    "35-39 years": ["35-39"],
    "40-44 years": ["40-44"],
    "45-49 years": ["45-49"],
    "50-54 years": ["50-54"],
    "55-59 years": ["55-59"],
    "60-64 years": ["60-64"],
    "65-69 years": ["65-69"],
    "70-74 years": ["70-74"],
    "75-79 years": ["75-79"],
    "80-84 years": ["80-84"],
    "85-89 years": ["85-89"],
    "90-94 years": ["90-94"],
    "95+ years": ["95-99", "100+"],
}


def wpp_pop(iso, sex_col, gbd_age, year):
    w = wfull[(wfull.Location == "World") if iso == "GLB" else (wfull.ISO3_code == iso)]
    w = w[(w.Time == year) & (w.AgeGrp.isin(AGE_MAP[gbd_age]))]
    return w[sex_col].sum()


pop_2053 = {}
all_iso = sorted(set(k[0] for k in pop_base.index))
for iso in all_iso:
    for sex, col in [("Male", "PopMale"), ("Female", "PopFemale")]:
        for ag in AGE15:
            p23, p53 = wpp_pop(iso, col, ag, 2023), wpp_pop(iso, col, ag, 2053)
            growth = p53 / p23 if p23 > 0 else 1.0
            pop_2053[(iso, sex, ag)] = float(pop_base.get((iso, sex, ag), 0)) * growth

our_iso3 = [i for i in all_iso if i != "GLB"]


def run_country(iso, pop_dict):
    total = {"diet": 0, "glp1": 0, "combo": 0}
    for o, fn in OUTFILES.items():
        inc = lt.load_measure(fn, "Incidence")
        prev_m = lt.load_measure(fn, "Prevalence")
        cd = lt.load_measure(fn, "Deaths")
        cfr_d, bg_d = {}, {}
        for k in set(inc.index) & set(prev_m.index) & set(cd.index):
            p = float(prev_m[k])
            cfr_d[k] = min(float(cd[k]) / p, 0.9) if p > 1e-6 else 0.0
            bg_d[k] = max(float(allc.get(k, 0)) - float(cd[k]), 0.0)
        dat = dict(inc=inc, prev=prev_m, cfr=cfr_d, bg=bg_d, pop=pop_dict)

        glp_row = f2[(f2.outcome == o) & (f2.iso3 == iso)]
        has_glp = len(glp_row) and glp_row.flow_available.iloc[0] in (True, "True")
        pif_med = float(glp_row.pif_med_pop.iloc[0]) if has_glp else 0.0
        pif_resid = float(glp_row.pif_resid_pop.iloc[0]) if has_glp else 0.0
        if o == "CKD":
            has_glp, pif_med, pif_resid = (
                True,
                0.0,
                0.24 * 0.20 * _ow_dict.get(iso, 0.5),
            )
        diet_pif = float(dp.get((o, iso), 0.0))
        RETAIN = 0.76

        pifs = {
            "baseline": (0.0, 0.0),
            "glp1": (pif_med * RETAIN, pif_resid * RETAIN),
            "diet": (diet_pif, 0.0),
            "combo": (diet_pif, pif_resid * RETAIN),
        }
        pifs2 = {"combo": pif_med * RETAIN}
        ramp2 = {"combo": 1}
        ramp_map = {
            "baseline": (99, 99),
            "glp1": (1, 1),
            "diet": (8, 8),
            "combo": (8, 1),
        }

        for sex in ["Male", "Female"]:
            sc = lt.simulate(iso, sex, o, dat, pifs, 30, ramp_map, pifs2, ramp2)
            for arm in ["glp1", "diet", "combo"]:
                if arm in sc and "baseline" in sc:
                    total[arm] += sum(b - v for b, v in zip(sc["baseline"], sc[arm]))
    return total


countries = our_iso3 + ["GLB"]
results = []
for iso in countries:
    print(f"Running {iso}...", flush=True)
    t23 = run_country(iso, pop_base)
    t53 = run_country(iso, pop_2053)

    d_chg = (t53["diet"] - t23["diet"]) / t23["diet"] * 100 if t23["diet"] > 0 else 0
    g_chg = (t53["glp1"] - t23["glp1"]) / t23["glp1"] * 100 if t23["glp1"] > 0 else 0
    r23 = t23["diet"] / t23["glp1"] if t23["glp1"] > 0 else np.inf
    r53 = t53["diet"] / t53["glp1"] if t53["glp1"] > 0 else np.inf
    r_chg = (r53 - r23) / r23 * 100 if r23 > 0 and r23 != np.inf else 0

    results.append(
        {
            "Country": iso,
            "Diet 2023 (k)": round(t23["diet"] / 1e3, 1),
            "Diet 2053 (k)": round(t53["diet"] / 1e3, 1),
            "Diet change (%)": round(d_chg, 1),
            "GLP-1 2023 (k)": round(t23["glp1"] / 1e3, 1),
            "GLP-1 2053 (k)": round(t53["glp1"] / 1e3, 1),
            "GLP-1 change (%)": round(g_chg, 1),
            "Ratio 2023": round(r23, 2),
            "Ratio 2053": round(r53, 2),
            "Ratio change (%)": round(r_chg, 1),
        }
    )

df = pd.DataFrame(results)
glb_row = df[df.Country == "GLB"]
rest = df[df.Country != "GLB"].sort_values("Country")
df = pd.concat([glb_row, rest])

with open("tables/tableS_population_sensitivity_v2.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Supplementary Table. Population sensitivity (2023 vs 2053)."])
    w.writerow(["Base case 2023 GBD population vs 2053 UN WPP medium variant."])
    w.writerow([])
    w.writerow(df.columns.tolist())
    for _, r in df.iterrows():
        w.writerow(r.tolist())

print(f"\n✓ Saved: {len(df)} rows")
glb = df[df.Country == "GLB"].iloc[0]
print(
    f"GLB: ratio {glb['Ratio 2023']}→{glb['Ratio 2053']} ({glb['Ratio change (%)']}%)"
)
print(f"Diet: {glb['Diet change (%)']}%, GLP-1: {glb['GLP-1 change (%)']}%")
print()
print(df.to_string(index=False))
