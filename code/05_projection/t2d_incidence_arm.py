import os
import pandas as pd, numpy as np, argparse

import pandas as _cmpd

_cm = _cmpd.read_csv("params/country_map.csv")
ISO = dict(zip(_cm.gbd_location_name, _cm.iso3))
AGE15 = [f"{a}-{a+4} years" for a in range(25, 95, 5)] + ["95+ years"]


def ramp(t, yrs):
    import math

    k = 4.0 / yrs
    return 1.0 / (1.0 + math.exp(-k * (t - yrs / 1.6)))


def main(a):
    d = pd.read_csv("data/gbd/cause_totals/queryB_t2d.csv", low_memory=False)
    d = d[
        (d.measure_name == "Incidence")
        & (d.metric_name == "Number")
        & (d.year == 2023)
        & d.age_name.isin(AGE15)
        & d.sex_name.isin(["Male", "Female"])
        & d.location_name.isin(ISO)
    ]
    d["iso3"] = d.location_name.map(ISO)
    inc0 = d.groupby("iso3").val.sum()

    tp = pd.read_csv("params/glp1_trial_params.csv")
    hr = float(tp[tp.param == "t2d_incidence_hr"].value.iloc[0])
    RETAIN_T2D = 0.72
    glp_pif = a.coverage * (1 - hr) * RETAIN_T2D
    _ow = pd.read_csv("params/overweight_prevalence.csv")
    ow = dict(zip(_ow.iso3, _ow.overweight_prev))

    dp = pd.read_csv(os.environ.get("LT_DP_OVERRIDE", "results/diet_policy_arm.csv"))
    dpq = dp[(dp.scenario == a.diet_scenario) & (dp.outcome == "T2D")]
    diet_pif = dpq.groupby("iso3").averted_frac_alldeaths.mean()

    rows = []
    for iso, base in inc0.items():
        for t in range(1, a.years + 1):
            rows.append(
                dict(
                    iso3=iso,
                    year=t,
                    glp1=base * glp_pif * ow.get(iso, 0.5) * ramp(t, a.ramp_glp1),
                    diet=base * float(diet_pif.get(iso, 0)) * ramp(t, a.ramp_diet),
                )
            )
    t2 = pd.DataFrame(rows)
    t2.to_csv(a.out, index=False)
    s = t2.groupby("iso3")[["glp1", "diet"]].sum() / 1e6
    print(f"=== New T2D cases averted, cumulative over {a.years} years (millions) ===")
    print(s.round(2).to_string())
    print(
        "\nNote: GLP-1 is an on-treatment effect (reverts on discontinuation); the dietary arm acts through sustained exposure change"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--coverage", type=float, default=0.20)
    p.add_argument("--diet-scenario", default="feasible")
    p.add_argument("--years", type=int, default=30)
    p.add_argument("--ramp-glp1", type=int, default=4)
    p.add_argument("--ramp-diet", type=int, default=8)
    p.add_argument("--out", default="results/t2d_incidence_arms.csv")
    main(p.parse_args())
