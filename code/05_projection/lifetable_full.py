import pandas as pd, numpy as np, argparse

import pandas as _cmpd

_cm = _cmpd.read_csv("params/country_map.csv")
ISO = dict(zip(_cm.gbd_location_name, _cm.iso3))
AGE15 = [f"{a}-{a+4} years" for a in range(25, 95, 5)] + ["95+ years"]
OUTFILES = {
    "IHD": "queryB_ihd.csv",
    "IschStroke": "queryB_istroke.csv",
    "ICH": "queryB_ich.csv",
    "SAH": "queryB_sah.csv",
    "CKD": "queryB_ckd_3years.csv",
}
MED_FRAC = 0.33


def load_measure(fname, measure, metric="Rate"):
    d = pd.read_csv(f"data/gbd/cause_totals/{fname}", low_memory=False)
    d = d[
        (d.measure_name == measure)
        & (d.metric_name == metric)
        & (d.year == 2023)
        & d.age_name.isin(AGE15)
        & d.sex_name.isin(["Male", "Female"])
        & d.location_name.isin(ISO)
    ]
    d["iso3"] = d.location_name.map(ISO)
    v = d.set_index(["iso3", "sex_name", "age_name"]).val
    return (v / 1e5) if metric == "Rate" else v


def load_pop():
    d = pd.read_csv("data/gbd/1-3_population.csv", low_memory=False)
    d = d[
        (d.year == 2023)
        & d.age_name.isin(AGE15)
        & d.sex_name.isin(["Male", "Female"])
        & d.location_name.isin(ISO)
    ]
    d["iso3"] = d.location_name.map(ISO)
    return d.set_index(["iso3", "sex_name", "age_name"]).val


def ramp(t, yrs):
    import math

    k = 4.0 / yrs
    return 1.0 / (1.0 + math.exp(-k * (t - yrs / 1.6)))


def simulate(iso, sex, o, dat, pifs, years, ramp_map, pifs2=None, ramp2=None):
    inc, prev, cfr, bg = (dat[k] for k in ["inc", "prev", "cfr", "bg"])
    g = lambda tbl, ag: float(tbl.get((iso, sex, ag), 0.0))
    scen = {}
    for name, (pif_inc, pif_cfr) in pifs.items():
        W = {ag: g(dat["pop"], ag) * (1 - g(prev, ag)) for ag in AGE15}
        D = {ag: g(dat["pop"], ag) * g(prev, ag) for ag in AGE15}
        deaths = []
        for t in range(1, years + 1):
            rt_inc = ramp(t, ramp_map[name][0])
            f2 = 1.0
            if pifs2 and name in pifs2:
                f2 = 1 - pifs2[name] * ramp(t, ramp2[name])
            rt_cfr = ramp(t, ramp_map[name][1]) * (0.5 if t == 1 else 1.0)
            yd = 0.0
            newW, newD = {}, {}
            for ag in AGE15:
                i_eff = (
                    g(inc, ag)
                    / max(1 - g(prev, ag), 1e-9)
                    * (1 - pif_inc * rt_inc)
                    * f2
                )
                f_eff = g(cfr, ag) * (1 - pif_cfr * rt_cfr)
                b = g(bg, ag)
                new_cases = W[ag] * i_eff
                d_cause = D[ag] * f_eff
                yd += d_cause
                newW[ag] = max(W[ag] - new_cases - W[ag] * b, 0.0)
                newD[ag] = max(D[ag] + new_cases - d_cause - D[ag] * b, 0.0)
            W2, D2 = {}, {}
            for i2, ag in enumerate(AGE15):
                W2[ag] = newW[ag] * 0.8 + (newW[AGE15[i2 - 1]] * 0.2 if i2 > 0 else 0)
                D2[ag] = newD[ag] * 0.8 + (newD[AGE15[i2 - 1]] * 0.2 if i2 > 0 else 0)
            W, D = W2, D2
            deaths.append(yd)
        scen[name] = deaths
    return scen


def main(a):
    pop = load_pop()
    allc = load_measure("queryD_All_causes_deaths.csv", "Deaths")
    import os as _os

    f2 = pd.read_csv(_os.environ.get("LT_FLOW_OVERRIDE", "results/figure2_flow_v2.csv"))
    f2 = f2[f2.coverage == 0.20]
    _ow = pd.read_csv("params/overweight_prevalence.csv")
    _ow_dict = dict(zip(_ow.iso3, _ow.overweight_prev))
    afford = (
        dict(pd.read_csv("params/coverage_afford_cap.csv").values)
        if a.coverage == "afford"
        else {}
    )
    dp = pd.read_csv(_os.environ.get("LT_DP_OVERRIDE", "results/diet_policy_arm.csv"))
    dp = (
        dp[dp.scenario == a.diet_scenario]
        .groupby(["outcome", "iso3"])
        .averted_frac_alldeaths.mean()
        / 100
        if dp[dp.scenario == a.diet_scenario].averted_frac_alldeaths.max() > 1
        else dp[dp.scenario == a.diet_scenario]
        .groupby(["outcome", "iso3"])
        .averted_frac_alldeaths.mean()
    )

    rows = []
    for o, fn in OUTFILES.items():
        inc = load_measure(fn, "Incidence")
        prev = load_measure(fn, "Prevalence")
        cd = load_measure(fn, "Deaths")
        cfr, bg = {}, {}
        for k in set(inc.index) & set(prev.index) & set(cd.index):
            p = float(prev[k])
            cfr[k] = min(float(cd[k]) / p, 0.9) if p > 1e-6 else 0.0
            bg[k] = max(float(allc.get(k, 0)) - float(cd[k]), 0.0)
        dat = dict(inc=inc, prev=prev, cfr=cfr, bg=bg, pop=pop)
        for iso in ISO.values():
            glp_key = (o, iso)
            glp_row = f2[(f2.outcome == o) & (f2.iso3 == iso)]
            has_glp = len(glp_row) and glp_row.flow_available.iloc[0] in (True, "True")
            pif_med = float(glp_row.pif_med_pop.iloc[0]) if has_glp else 0.0
            pif_resid = float(glp_row.pif_resid_pop.iloc[0]) if has_glp else 0.0
            _split = _os.environ.get("SENS_SPLIT_INC")
            if _split is not None and o in ("IHD", "IschStroke") and has_glp:
                _tot = pif_med + pif_resid
                pif_med, pif_resid = float(_split) * _tot, (1 - float(_split)) * _tot
            _hr = _os.environ.get("SENS_CV_HR")
            if _hr is not None and o in ("IHD", "IschStroke") and has_glp:
                _k = (1 - float(_hr)) / 0.20
                pif_med, pif_resid = pif_med * _k, pif_resid * _k
            if o == "CKD":
                has_glp, pif_med, pif_resid = (
                    True,
                    0.0,
                    (1 - float(_os.environ.get("SENS_CKD_HR", 0.76)))
                    * 0.20
                    * _ow_dict.get(iso, 0.5),
                )
            _ef = _os.environ.get("SENS_ELIG_FILE")
            if _ef:
                _sc = dict(pd.read_csv(_ef).values).get(iso, 1.0)
                pif_med, pif_resid = pif_med * _sc, pif_resid * _sc
            cov_eff = (
                afford.get(iso, 1.0) if a.coverage == "afford" else float(a.coverage)
            )
            pif_med *= cov_eff / 0.20
            pif_resid *= cov_eff / 0.20
            diet_pif = float(dp.get((o, iso), 0.0))
            import os

            RETAIN = float(os.environ.get("SENS_RETAIN", 0.76))
            RETAIN_T2D = float(os.environ.get("SENS_RETAIN_T2D", 0.72))
            if o == "T2D":
                RETAIN = RETAIN_T2D
            PM = float(os.environ.get("SENS_POLICY_MAINTAIN", 1.0))
            pifs = {
                "baseline": (0.0, 0.0),
                "glp1": (pif_med * RETAIN, pif_resid * RETAIN),
                "diet": (diet_pif * PM, 0.0),
                "combo": (diet_pif * PM, pif_resid * RETAIN),
            }
            pifs2 = {"combo": pif_med * RETAIN}
            ramp2 = {"combo": a.ramp_glp1}
            ramp_map = {
                "baseline": (99, 99),
                "glp1": (a.ramp_glp1, a.ramp_glp1),
                "diet": (a.ramp_diet, a.ramp_diet),
                "combo": (a.ramp_diet, a.ramp_glp1),
            }
            for sex in ["Male", "Female"]:
                sc = simulate(iso, sex, o, dat, pifs, a.years, ramp_map, pifs2, ramp2)
                for name in ["glp1", "diet", "combo"]:
                    if name == "glp1" and not has_glp:
                        continue
                    if name == "diet" and diet_pif == 0.0:
                        continue
                    if name == "combo" and not (has_glp or diet_pif > 0):
                        continue
                    averted = [b - v for b, v in zip(sc["baseline"], sc[name])]
                    row = dict(
                        iso3=iso,
                        sex=sex,
                        outcome=o,
                        arm=name,
                        averted_10yr=sum(averted),
                        averted_yr1=averted[0],
                        averted_yr10=averted[-1],
                    )
                    for yi, v in enumerate(averted, 1):
                        row[f"yr{yi}"] = v
                    rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(a.out, index=False)
    print(
        f"=== Full method: cumulative avertable deaths over {a.years} years (thousands, by arm) ==="
    )
    print(
        (t.groupby(["arm", "iso3"]).averted_10yr.sum() / 1000)
        .round(1)
        .unstack(0)
        .to_string()
    )
    print("\n=== Lag shape: year 1 vs year 10 (IHD, whole population) ===")
    q = (
        t[t.outcome == "IHD"]
        .groupby(["arm", "iso3"])[["averted_yr1", "averted_yr10"]]
        .sum()
    )
    print(q.round(0).to_string())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--years", type=int, default=10)
    p.add_argument("--ramp-glp1", type=int, default=1)
    p.add_argument("--ramp-diet", type=int, default=8)
    p.add_argument("--diet-scenario", default="feasible")
    p.add_argument("--coverage", default="0.20")
    p.add_argument("--out", default="results/figure4_full.csv")
    main(p.parse_args())
