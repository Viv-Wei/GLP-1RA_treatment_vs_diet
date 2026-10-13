import os
import pandas as pd, numpy as np, glob

COV = [0.05, 0.10, 0.20, 0.50]
import pandas as _cmpd

_cm = _cmpd.read_csv("params/country_map.csv")
ISO = dict(zip(_cm.gbd_location_name, _cm.iso3))
AGE15 = [f"{a}-{a+4} years" for a in range(25, 95, 5)] + ["95+ years"]
MED_FRAC = 0.33


def official_paf_diet(glob_pat):
    g = glob.glob(glob_pat)
    d = pd.read_csv(g[0], low_memory=False)
    d = d[
        (d.year == 2023)
        & (d.sex_name == "Both")
        & d.age_name.isin(AGE15)
        & d.location_name.isin(ISO)
        & (d.measure_name == "Deaths")
    ]
    num = d[d.metric_name == "Number"].set_index(["location_name", "age_name"]).val
    pct = d[d.metric_name == "Percent"].set_index(["location_name", "age_name"]).val
    attr = num.groupby("location_name").sum()
    tot = (num / pct.replace(0, np.nan)).groupby("location_name").sum()
    return (attr / tot).rename(index=ISO)


def merged_paf_diet(path):
    d = pd.read_csv(path, low_memory=False)
    d = d[
        (d.rei_name == "Dietary risks")
        & (d.year == 2023)
        & (d.sex_name == "Both")
        & d.age_name.isin(AGE15)
        & d.location_name.isin(ISO)
        & (d.measure_name == "Deaths")
    ]
    num = d[d.metric_name == "Number"].set_index(["location_name", "age_name"]).val
    pct = d[d.metric_name == "Percent"].set_index(["location_name", "age_name"]).val
    attr = num.groupby("location_name").sum()
    tot = (num / pct.replace(0, np.nan)).groupby("location_name").sum()
    return (attr / tot).rename(index=ISO)


PAF = {
    "IHD": official_paf_diet("data/gbd/Ischemic heart disease/Dietary risks/*/*.csv"),
    "IschStroke": official_paf_diet("data/gbd/Ischemic stroke/Dietary risks/*/*.csv"),
    "T2D": official_paf_diet("data/gbd/Diabetes mellitus type 2/Dietary risks/*/*.csv"),
    "ICH": merged_paf_diet("data/gbd/queryA_ICH_all_merged.csv"),
    "SAH": merged_paf_diet("data/gbd/queryA_SAH_all_merged.csv"),
}

f1 = pd.read_csv(os.environ.get("CRA_F1CH", "results/figure1_channel_burden.csv"))


def pool(outcome):
    w = (
        f1[f1.outcome == outcome]
        .pivot_table(
            index="iso3", columns="glp1_channel", values="deaths", aggfunc="sum"
        )
        .fillna(0)
    )
    return w.get("Full convergence - structurally reachable", 0) / w.clip(lower=0).sum(
        axis=1
    )


pif = pd.read_csv(os.environ.get("CRA_GLP1PIF", "results/glp1_arm_pif.csv"))
HR_TOT = 0.80
PIFCOL = {
    "IHD": ["pif_ihd_sbp", "pif_ihd_fpg"],
    "IschStroke": ["pif_str_sbp", "pif_is_fpg"],
}

_ow = pd.read_csv("params/overweight_prevalence.csv")
_ow_dict = dict(zip(_ow.iso3, _ow.overweight_prev))

rows = []
for outc in ["IHD", "IschStroke", "T2D", "ICH", "SAH"]:
    sr = pool(outc)
    for iso in sr.index:
        for c in COV:
            r = dict(
                outcome=outc,
                iso3=iso,
                coverage=c,
                PAF_diet=float(PAF[outc].get(iso, np.nan)),
                share_reach=float(sr[iso]),
                flow_available=outc in PIFCOL or outc == "T2D",
            )
            if outc in PIFCOL:
                r["pif_med_pop"] = c * _ow_dict.get(iso, 0.5) * MED_FRAC * (1 - HR_TOT)
                r["pif_resid_pop"] = (
                    c * _ow_dict.get(iso, 0.5) * (1 - MED_FRAC) * (1 - HR_TOT)
                )
            elif outc == "T2D":
                r["pif_med_pop"] = np.nan
                r["pif_resid_pop"] = np.nan
            else:
                r["pif_med_pop"] = np.nan
                r["pif_resid_pop"] = np.nan
            rows.append(r)

t = pd.DataFrame(rows)
t["averted_med_raw"] = t.pif_med_pop / t.PAF_diet
t["averted_med_samebase"] = t[["averted_med_raw", "share_reach"]].min(axis=1)
t["spillover_nondiet"] = (t.averted_med_raw - t.averted_med_samebase).clip(lower=0)
t["averted_resid_samebase"] = t.pif_resid_pop / t.PAF_diet
t["flow_ratio"] = (t.averted_med_samebase / t.share_reach).clip(upper=1.0)
t.to_csv(os.environ.get("CRA_F2FLOW_OUT", "results/figure2_flow_v2.csv"), index=False)

v = t[t.coverage == 0.20]
print("=== 20% coverage: five outcomes x six regions ===")
show = v.pivot_table(
    index="iso3",
    columns="outcome",
    values=["share_reach", "averted_med_samebase", "flow_ratio"],
)
for m, lab in [
    ("share_reach", "reachable pool %"),
    ("averted_med_samebase", "mediated leverage %"),
    ("flow_ratio", "flow ratio"),
]:
    print(f"\n--- {lab} ---")
    q = show[m][["IHD", "IschStroke", "ICH", "SAH", "T2D"]] * (
        100 if m != "flow_ratio" else 1
    )
    print(q.round(2 if m == "flow_ratio" else 1).to_string())
print("\n-> results/figure2_flow_v2.csv")
