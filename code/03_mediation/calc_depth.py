import os, json, argparse
import pandas as pd, numpy as np
from importlib.machinery import SourceFileLoader

pv2 = SourceFileLoader("pv2", "code/03_mediation/paf_deficiency_v2.py").load_module()

_cm = pd.read_csv("params/country_map.csv")
NAME2ISO = dict(zip(_cm.gbd_location_name, _cm.iso3))
_bl = pd.read_csv("params/mediator_baselines.csv")


def baseline(iso3, med):
    b = _bl[(_bl.iso3 == iso3) & (_bl.mediator == med)]
    both = b[b.sex == "Both"]
    return float(both.baseline.iloc[0]) if len(both) else float(b.baseline.mean())


RISK2FACTOR = {
    "Diet low in whole grains": "Whole grains (low)",
    "Diet low in fruits": "Fruits (low)",
    "Diet low in vegetables": "Vegetables (low)",
    "Diet low in fibre": "Fibre (low)",
    "Diet low in legumes": "Legumes (low)",
    "Diet low in nuts and seeds": "Nuts & seeds (low)",
    "Diet low in seafood omega-3 fatty acids": "Seafood omega-3 (low)",
    "Diet low in omega-6 polyunsaturated fatty acids": "Omega-6 PUFA (low)",
    "Diet high in red meat": "Red meat (high)",
    "Diet high in processed meat": "Processed meat (high)",
    "Diet high in sugar-sweetened beverages": "SSB (high)",
    "Diet high in sodium": "Sodium (high)",
    "Diet high in trans fatty acids": "Trans fat (high)",
}
MED2KEY = {
    "High systolic blood pressure": "SBP",
    "High fasting plasma glucose": "FPG",
    "High LDL cholesterol": "LDL",
}
SBP_CURVE = {
    "IHD": "high_sbp_ihd",
    "IschStroke": "high_sbp_s",
    "ICH": "high_sbp_s",
    "SAH": "high_sbp_s",
    "CKD": "high_sbp_ckd",
}
FPG_CURVE = {
    "IHD": "high_glucose_ihd",
    "IschStroke": "high_glucose_istroke",
    "ICH": "high_glucose_istroke",
}
TOTFILE = {
    "IHD": "queryB_ihd.csv",
    "IschStroke": "queryB_istroke.csv",
    "ICH": "queryB_ich.csv",
    "SAH": "queryB_sah.csv",
    "CKD": "queryB_ckd_3years.csv",
    "T2D": "queryB_t2d.csv",
}
B_TRIAL = {"T2D": {"FPG": 0.146}, "CKD": {"FPG": 0.048}}
COV_BASE = 0.20
_ow = pd.read_csv("params/overweight_prevalence.csv")
_ow_dict = dict(zip(_ow.iso3, _ow.overweight_prev))
GLP_SBP, GLP_FPG, GLP_LDL_REL = 5.10, 0.437, 0.04
LDL_RR_PER_MMOL = 1 / 0.78


def shift_from_rr(cur, base, rr, span, n=1000):
    for sp in (span, span * 0.75, span * 0.5, span * 0.25):
        try:
            grid = np.linspace(0, sp, n)
            rrs = cur.rr(base + grid) / cur.rr(base)
            break
        except ValueError:
            continue
    else:
        raise
    if rr > rrs[-1]:
        print(
            f"[warning] inversion: RR {rr:.3f} exceeds grid upper bound {rrs[-1]:.3f} (span={sp}), using upper bound"
        )
    return float(np.interp(rr, rrs, grid))


def total_deaths(iso3, outcome):
    d = pd.read_csv(f"data/gbd/cause_totals/{TOTFILE[outcome]}", low_memory=False)
    d = d[
        (d.measure_name == "Deaths")
        & (d.metric_name == "Number")
        & (d.year == 2023)
        & d.sex_name.isin(["Male", "Female"])
        & d.age_name.str.contains("years")
        & ~d.age_name.str.contains("standard")
    ]
    d = d[d.location_name.map(NAME2ISO) == iso3]
    return float(d.val.sum())


def main(a):
    dr = pd.read_csv(
        os.environ.get("CRA_CURVES", "params/dose_response_gbd2023_gbdexp.csv")
    )
    md = pd.read_csv("params/mediation_gbd2023_5outcomes.csv")
    md = md[md.outcome == a.outcome]
    flows = pd.read_csv(f"results/figure3_flows_{a.iso3}_{a.outcome}.csv")
    fac_deaths = flows.groupby("factor").deaths.sum()
    tot = total_deaths(a.iso3, a.outcome)

    sbp_cur = (
        pv2.Curve(dr[dr.curve_id == SBP_CURVE[a.outcome]])
        if a.outcome in SBP_CURVE
        else None
    )
    fpg_cur = (
        pv2.Curve(dr[dr.curve_id == FPG_CURVE[a.outcome]])
        if a.outcome in FPG_CURVE
        else None
    )
    sbp_base = baseline(a.iso3, "SBP")

    push = {"SBP": 0.0, "FPG": 0.0, "LDL": 0.0}
    for _, r in md.iterrows():
        m = MED2KEY.get(r.mediator)
        fac = RISK2FACTOR.get(r.risk)
        if not m or not fac:
            continue
        paf = fac_deaths.get(fac, 0.0) / tot
        if paf <= 0:
            continue
        rr_f = (1.0 / (1.0 - paf)) ** float(r.mf)
        if m == "SBP" and sbp_cur is not None:
            push["SBP"] += shift_from_rr(sbp_cur, sbp_base, rr_f, 40)
        elif m == "FPG" and fpg_cur is not None:
            push["FPG"] += shift_from_rr(fpg_cur, 5.1, rr_f, 4)
        elif m == "LDL":
            push["LDL"] += float(np.log(rr_f) / np.log(LDL_RR_PER_MMOL))

    COV = COV_BASE * _ow_dict.get(a.iso3, 0.5)
    GLP = {"SBP": GLP_SBP, "FPG": GLP_FPG, "LDL": GLP_LDL_REL * baseline(a.iso3, "LDL")}
    depth = {}
    for m in GLP:
        if a.outcome in B_TRIAL and m == "FPG":
            depth[m] = B_TRIAL[a.outcome]["FPG"]
        elif push[m] <= 0:
            depth[m] = -1.0
        else:
            depth[m] = round(min(COV * GLP[m] / push[m], 1.0), 3)
    os.makedirs("params", exist_ok=True)
    json.dump(
        {
            **depth,
            "push": {k: round(v, 4) for k, v in push.items()},
            "glp_shift": {k: round(COV * v, 4) for k, v in GLP.items()},
        },
        open(f"params/figure3_depth_{a.iso3}_{a.outcome}.json", "w"),
    )
    print(
        a.iso3,
        a.outcome,
        "total deaths:",
        round(tot),
        "inflow:",
        {k: round(v, 3) for k, v in push.items()},
        "depth:",
        depth,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--iso3")
    p.add_argument("--outcome")
    main(p.parse_args())
