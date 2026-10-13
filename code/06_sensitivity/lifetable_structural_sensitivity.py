import pandas as pd, numpy as np, csv, os, sys
from multiprocessing import Pool

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))
from importlib.machinery import SourceFileLoader

lt = SourceFileLoader("lt", "code/05_projection/lifetable_full.py").load_module()

pop = lt.load_pop()
allc = lt.load_measure("queryD_All_causes_deaths.csv", "Deaths")
f2 = pd.read_csv(os.environ.get("LT_FLOW_OVERRIDE", "results/figure2_flow_v2.csv"))
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
RETAIN = 0.76
OUTFILES = {
    "IHD": "queryB_ihd.csv",
    "IschStroke": "queryB_istroke.csv",
    "CKD": "queryB_ckd_3years.csv",
}

DAT = {}
for o, fn in OUTFILES.items():
    inc = lt.load_measure(fn, "Incidence")
    prev_m = lt.load_measure(fn, "Prevalence")
    cd = lt.load_measure(fn, "Deaths")
    cfr_d, bg_d = {}, {}
    for k in set(inc.index) & set(prev_m.index) & set(cd.index):
        p = float(prev_m[k])
        cfr_d[k] = min(float(cd[k]) / p, 0.9) if p > 1e-6 else 0.0
        bg_d[k] = max(float(allc.get(k, 0)) - float(cd[k]), 0.0)
    DAT[o] = dict(inc=inc, prev=prev_m, cfr=cfr_d, bg=bg_d, pop=pop)


def simulate_struct(
    iso,
    sex,
    dat,
    pifs,
    years,
    ramp_map,
    pifs2,
    ramp2,
    half_cycle_frac=0.0,
    bg_mult_D=1.0,
):
    inc, prev, cfr, bg = (dat[k] for k in ["inc", "prev", "cfr", "bg"])
    g = lambda tbl, ag: float(tbl.get((iso, sex, ag), 0.0))
    scen = {}
    for name, (pif_inc, pif_cfr) in pifs.items():
        W = {ag: g(dat["pop"], ag) * (1 - g(prev, ag)) for ag in lt.AGE15}
        D = {ag: g(dat["pop"], ag) * g(prev, ag) for ag in lt.AGE15}
        deaths = []
        for t in range(1, years + 1):
            rt_inc = lt.ramp(t, ramp_map[name][0])
            f2v = (
                1 - pifs2[name] * lt.ramp(t, ramp2[name])
                if (pifs2 and name in pifs2)
                else 1.0
            )
            rt_cfr = lt.ramp(t, ramp_map[name][1]) * (0.5 if t == 1 else 1.0)
            yd, newW, newD = 0.0, {}, {}
            for ag in lt.AGE15:
                i_eff = (
                    g(inc, ag)
                    / max(1 - g(prev, ag), 1e-9)
                    * (1 - pif_inc * rt_inc)
                    * f2v
                )
                f_eff = g(cfr, ag) * (1 - pif_cfr * rt_cfr)
                b_W = g(bg, ag)
                b_D = b_W * bg_mult_D
                new_cases = W[ag] * i_eff
                d_cause = D[ag] * f_eff + half_cycle_frac * new_cases * f_eff
                yd += d_cause
                newW[ag] = max(W[ag] - new_cases - W[ag] * b_W, 0.0)
                newD[ag] = max(D[ag] + new_cases - d_cause - D[ag] * b_D, 0.0)
            W2, D2 = {}, {}
            for i2, ag in enumerate(lt.AGE15):
                W2[ag] = newW[ag] * 0.8 + (
                    newW[lt.AGE15[i2 - 1]] * 0.2 if i2 > 0 else 0
                )
                D2[ag] = newD[ag] * 0.8 + (
                    newD[lt.AGE15[i2 - 1]] * 0.2 if i2 > 0 else 0
                )
            W, D = W2, D2
            deaths.append(yd)
        scen[name] = deaths
    return scen


def run(args):
    iso, scen_name, hcf, bgm = args
    total = {"diet": 0.0, "glp1": 0.0, "combo": 0.0}
    for o in OUTFILES:
        glp_row = f2[(f2.outcome == o) & (f2.iso3 == iso) & (f2.coverage == 0.20)]
        has_glp = len(glp_row) and glp_row.flow_available.iloc[0] in (True, "True")
        pif_med = float(glp_row.pif_med_pop.iloc[0]) if has_glp else 0.0
        pif_resid = float(glp_row.pif_resid_pop.iloc[0]) if has_glp else 0.0
        if o == "CKD":
            pif_med, pif_resid = 0.0, 0.24 * 0.20 * _ow_dict.get(iso, 0.5)
        diet_pif = float(dp.get((o, iso), 0.0))
        pifs = {
            "baseline": (0.0, 0.0),
            "glp1": (pif_med * RETAIN, pif_resid * RETAIN),
            "diet": (diet_pif, 0.0),
            "combo": (diet_pif, pif_resid * RETAIN),
        }
        pifs2, ramp2 = {"combo": pif_med * RETAIN}, {"combo": 1}
        ramp_map = {
            "baseline": (99, 99),
            "glp1": (1, 1),
            "diet": (8, 8),
            "combo": (8, 1),
        }
        for sex in ["Male", "Female"]:
            sc = simulate_struct(
                iso, sex, DAT[o], pifs, 30, ramp_map, pifs2, ramp2, hcf, bgm
            )
            for arm in total:
                total[arm] += sum(b - v for b, v in zip(sc["baseline"], sc[arm]))
    return dict(iso3=iso, scenario=scen_name, **total)


SCEN = [
    ("Base case", 0.0, 1.0),
    ("Same-cycle CFR for 50% of new cases", 0.5, 1.0),
    ("Background mortality 1.5x for diseased", 0.0, 1.5),
    ("Background mortality 2.0x for diseased", 0.0, 2.0),
    ("Background mortality 3.0x for diseased", 0.0, 3.0),
]

if __name__ == "__main__":
    isos = sorted(set(lt.ISO.values()) | {"GLB"})
    jobs = [(iso, n, h, b) for iso in isos for n, h, b in SCEN]
    with Pool(min(len(jobs), 100)) as p:
        res = pd.DataFrame(p.map(run, jobs))
    res["ratio"] = res.diet / res.glp1
    os.makedirs("results/itt_1012", exist_ok=True)
    res.to_csv("results/itt_1012/sens_structural_countries_itt.csv", index=False)

    g = res[res.iso3 == "GLB"].set_index("scenario").loc[[s[0] for s in SCEN]]
    b = g.loc["Base case"]
    df = pd.DataFrame(
        {
            "Scenario": g.index,
            "Diet (M)": (g.diet / 1e6).round(2).values,
            "GLP-1 (M)": (g.glp1 / 1e6).round(2).values,
            "Combo (M)": (g.combo / 1e6).round(2).values,
            "Diet/GLP-1": g.ratio.round(2).values,
            "Diet change (%)": ((g.diet / b.diet - 1) * 100).round(1).values,
            "GLP-1 change (%)": ((g.glp1 / b.glp1 - 1) * 100).round(1).values,
            "Ratio change (%)": ((g.ratio / b.ratio - 1) * 100).round(1).values,
        }
    )
    with open(
        "tables/tableS_lifetable_structural_sensitivity.csv", "w", newline=""
    ) as f:
        w = csv.writer(f)
        w.writerow(
            [
                "Supplementary Table. Sensitivity of life table results to structural assumptions (Global)."
            ]
        )
        w.writerow(
            [
                "Same-cycle case fatality: 50% of new incident cases exposed to CFR in cycle of onset. "
                "Background mortality: diseased compartment background mortality multiplied by 1.5, 2.0, or 3.0. "
                "All scenarios at 20% coverage with the feasible dietary package, Global."
            ]
        )
        w.writerow([])
        w.writerow(df.columns.tolist())
        for _, row in df.iterrows():
            w.writerow(row.tolist())
    print(df.to_string(index=False))
    lead = (
        res[res.iso3 != "GLB"]
        .groupby("scenario")
        .ratio.apply(lambda x: int((x > 1).sum()))
    )
    print("\nNumber of countries where diet leads:\n" + lead.to_string())
    print("\nBase case should be about diet 9.96M, glp1 4.06M, ratio 2.45")
