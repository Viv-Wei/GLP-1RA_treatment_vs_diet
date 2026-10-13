import os


import pandas as pd, numpy as np
from importlib.machinery import SourceFileLoader

dpa = SourceFileLoader(
    "dpa", "code/02_intervention_arms/diet_policy_arm.py"
).load_module()
pv2 = dpa.pv2

import sys

SHARD = int(sys.argv[sys.argv.index("--shard") + 1]) if "--shard" in sys.argv else 0
NSH = int(sys.argv[sys.argv.index("--nshards") + 1]) if "--nshards" in sys.argv else 1
RNG = np.random.default_rng(42 + SHARD)
N = 1000 // NSH
MED_FRAC = 0.33
_ow = pd.read_csv("params/overweight_prevalence.csv")
_ow_dict = dict(zip(_ow.iso3, _ow.overweight_prev))


def draw_normal(mean, lo, hi, n):
    return RNG.normal(mean, (hi - lo) / (2 * 1.96), n)


tp = pd.read_csv("params/glp1_trial_params.csv")


def par(param, trial):
    r = tp[(tp.param == param) & (tp.trial == trial)].iloc[0]
    return float(r.value), float(r.ci_lo), float(r.ci_hi)


hr_m, hr_lo, hr_hi = par("mace_hr", "SELECT")
HR = draw_normal(np.log(hr_m), np.log(hr_lo), np.log(hr_hi), N)
HR = np.exp(HR)
NA2SBP = draw_normal(3.82 / 2.30, 3.08 / 2.30, 4.55 / 2.30, N)
AMP = RNG.uniform(1.5, 2.2, N)

curves = pd.read_csv(
    os.environ.get("CRA_CURVES", "params/dose_response_gbd2023_gbdexp.csv")
)
cvt = pd.read_csv("params/cv_by_factor.csv")
paf_f2 = pd.read_csv(os.environ.get("CRA_F2FLOW", "results/figure2_flow_v2.csv"))
paf_diet = (
    paf_f2[(paf_f2.outcome.isin(["IHD", "IschStroke"])) & (paf_f2.coverage == 0.20)]
    .set_index(["outcome", "iso3"])
    .PAF_diet
)

expo = {s: pd.read_csv(f"params/exposure_gbd_{s}.csv") for s in ["male", "female"]}
sbpf = {s: f"params/sbp_gbd2023_{s}.csv" for s in ["male", "female"]}

res = {
    k: []
    for k in ["sodium30", "feasible", "aspirational", "glp1_total", "glp1_mediated"]
}
for i in range(N):
    dpa.NA2SBP = NA2SBP[i]
    cvmap = dict(zip(cvt.gdd_var, cvt.cv_between * AMP[i]))
    per = {}
    SBP_CURVE = {"IHD": "high_sbp_ihd", "IschStroke": "high_sbp_s"}
    for s in ["male", "female"]:
        e = expo[s]
        base = dpa.combined_paf(e, curves, cvmap, outcomes=("IHD", "IschStroke"))
        nb = {
            o: dpa.sodium_paf_via_sbp(e, curves, sbpf[s], curve_id=c)
            for o, c in SBP_CURVE.items()
        }
        for scen in ["sodium30", "feasible", "aspirational"]:
            sh = dpa.shifted_exposure(e, curves, scen)
            sc = dpa.combined_paf(sh, curves, cvmap, outcomes=("IHD", "IschStroke"))
            ns = {
                o: dpa.sodium_paf_via_sbp(sh, curves, sbpf[s], curve_id=c)
                for o, c in SBP_CURVE.items()
            }
            for outc, iso in base:
                b = 1 - (1 - base[(outc, iso)]) * (1 - nb[outc][iso])
                v = 1 - (1 - sc[(outc, iso)]) * (1 - ns[outc][iso])
                per.setdefault((scen, outc, iso), []).append(b - v)
    for (scen, outc, iso), vals in per.items():
        res.setdefault(scen, []).append(
            (
                SHARD * N + i,
                outc,
                iso,
                np.mean(vals) * dpa.calib(outc, iso) / paf_diet[(outc, iso)] * 100,
            )
        )
    for outc, iso in paf_diet.index:
        _ow_iso = _ow_dict.get(iso, 0.5)
        pif_t = 0.20 * _ow_iso * (1 - HR[i])
        pif_m = 0.20 * _ow_iso * MED_FRAC * (1 - HR[i])
        res["glp1_total"].append(
            (SHARD * N + i, outc, iso, pif_t / paf_diet[(outc, iso)] * 100)
        )
        res["glp1_mediated"].append(
            (SHARD * N + i, outc, iso, pif_m / paf_diet[(outc, iso)] * 100)
        )
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{N}")

rows = []
draws = []
for scen, lst in res.items():
    df = pd.DataFrame(lst, columns=["draw", "outcome", "iso3", "v"])
    df["scenario"] = scen
    draws.append(df)
    q = df.groupby(["outcome", "iso3"]).v.quantile([0.5, 0.025, 0.975]).unstack()
    q.columns = ["med", "lo", "hi"]
    q["scenario"] = scen
    rows.append(q.reset_index())
out = pd.concat(rows)
out.to_csv(f"results/figure2b_arms_ui_19c_shard{SHARD}.csv", index=False)
pd.concat(draws).to_csv(f"results/f2b_draws_shard{SHARD}.csv", index=False)
print(out.round(1).to_string(index=False))
