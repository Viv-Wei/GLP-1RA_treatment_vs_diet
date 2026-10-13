import os, argparse
import pandas as pd, numpy as np
from multiprocessing import Pool
from importlib.machinery import SourceFileLoader

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))

ap = argparse.ArgumentParser()
ap.add_argument("--workers", type=int, default=100)
ap.add_argument("--n-ckd", type=int, default=200)
ap.add_argument("--n-t2d", type=int, default=100)
ap.add_argument("--suffix", default="")
A = ap.parse_args()

_ow = pd.read_csv("params/overweight_prevalence.csv")
_ow_dict = dict(zip(_ow.iso3, _ow.overweight_prev))
dpa = SourceFileLoader(
    "dpa", "code/02_intervention_arms/diet_policy_arm.py"
).load_module()

RNG = np.random.default_rng(99)


def draw_normal(mean, lo, hi, n):
    return RNG.normal(mean, (hi - lo) / (2 * 1.96), n)


tp = pd.read_csv("params/glp1_trial_params.csv")


def par(param, trial):
    r = tp[(tp.param == param) & (tp.trial == trial)].iloc[0]
    return float(r.value), float(r.ci_lo), float(r.ci_hi)


N_MAX = 200
hr_m, hr_lo, hr_hi = par("mace_hr", "SELECT")
HR_SELECT = np.exp(draw_normal(np.log(hr_m), np.log(hr_lo), np.log(hr_hi), N_MAX))
HR_CKD = np.exp(draw_normal(np.log(0.76), np.log(0.66), np.log(0.88), N_MAX))
t2d_m, t2d_lo, t2d_hi = par("t2d_incidence_hr", "SELECT_glycemia")
HR_T2D = np.exp(draw_normal(np.log(t2d_m), np.log(t2d_lo), np.log(t2d_hi), N_MAX))
NA2SBP_draws = draw_normal(3.82 / 2.30, 3.08 / 2.30, 4.55 / 2.30, N_MAX)
AMP_draws = RNG.uniform(1.5, 2.2, N_MAX)

curves = pd.read_csv(
    os.environ.get("CRA_CURVES", "params/dose_response_gbd2023_gbdexp.csv")
)
cvt = pd.read_csv("params/cv_by_factor.csv")
paf_f2 = pd.read_csv(os.environ.get("CRA_F2FLOW", "results/figure2_flow_v2.csv"))
paf_ckd = (
    paf_f2[(paf_f2.outcome == "CKD") & (paf_f2.coverage == 0.20)]
    .set_index(["outcome", "iso3"])
    .PAF_diet
)
if paf_ckd.empty:
    _dp = pd.read_csv(os.environ.get("LT_DP_OVERRIDE", "results/diet_policy_arm.csv"))
    _dp = _dp[(_dp.outcome == "CKD") & (_dp.scenario == "feasible")]
    _p = _dp.groupby("iso3").paf_base.mean()
    paf_ckd = pd.Series(
        _p.values,
        index=pd.MultiIndex.from_product(
            [["CKD"], _p.index], names=["outcome", "iso3"]
        ),
    )
    print(
        f"CKD diet PAF (sodium pathway) from dietary policy table: {len(paf_ckd)} locations",
        flush=True,
    )
paf_t2d = (
    paf_f2[(paf_f2.outcome == "T2D") & (paf_f2.coverage == 0.20)]
    .set_index(["outcome", "iso3"])
    .PAF_diet
)
expo = {s: pd.read_csv(f"params/exposure_gbd_{s}.csv") for s in ["male", "female"]}
sbpf = {s: f"params/sbp_gbd2023_{s}.csv" for s in ["male", "female"]}
SCENS = ["sodium30", "feasible", "aspirational"]


def ckd_draw(i):
    dpa.NA2SBP = NA2SBP_draws[i]
    rows, per = [], {}
    for s in ["male", "female"]:
        e = expo[s]
        nb_base = dpa.sodium_paf_via_sbp(e, curves, sbpf[s], curve_id="high_sbp_ckd")
        for scen in SCENS:
            sh = dpa.shifted_exposure(e, curves, scen)
            nb_scen = dpa.sodium_paf_via_sbp(
                sh, curves, sbpf[s], curve_id="high_sbp_ckd"
            )
            for iso in nb_base:
                per.setdefault((scen, iso), []).append(nb_base[iso] - nb_scen[iso])
    for (scen, iso), vals in per.items():
        if ("CKD", iso) in paf_ckd.index and paf_ckd[("CKD", iso)] > 0:
            rows.append(
                (
                    i,
                    "CKD",
                    iso,
                    scen,
                    np.mean(vals) * dpa.calib("CKD", iso) / paf_ckd[("CKD", iso)] * 100,
                )
            )
    pif = 0.20 * (1 - HR_CKD[i])
    for iso in paf_ckd.index.get_level_values("iso3"):
        if paf_ckd[("CKD", iso)] > 0:
            rows.append(
                (
                    i,
                    "CKD",
                    iso,
                    "glp1_total",
                    pif * _ow_dict.get(iso, 0.5) / paf_ckd[("CKD", iso)] * 100,
                )
            )
    return rows


def t2d_draw(i):
    dpa.NA2SBP = NA2SBP_draws[i]
    cvmap = dict(zip(cvt.gdd_var, cvt.cv_between * AMP_draws[i]))
    rows, per = [], {}
    for s in ["male", "female"]:
        e = expo[s]
        base = dpa.combined_paf(e, curves, cvmap, outcomes=("T2D",))
        for scen in SCENS:
            sh = dpa.shifted_exposure(e, curves, scen)
            sc = dpa.combined_paf(sh, curves, cvmap, outcomes=("T2D",))
            for outc, iso in base:
                per.setdefault((scen, iso), []).append(
                    base[(outc, iso)] - sc.get((outc, iso), 0.0)
                )
    for (scen, iso), vals in per.items():
        if ("T2D", iso) in paf_t2d.index and paf_t2d[("T2D", iso)] > 0:
            rows.append(
                (
                    i,
                    "T2D",
                    iso,
                    scen,
                    np.mean(vals) * dpa.calib("T2D", iso) / paf_t2d[("T2D", iso)] * 100,
                )
            )
    pif = 0.20 * (1 - HR_T2D[i])
    for iso in paf_t2d.index.get_level_values("iso3"):
        if paf_t2d[("T2D", iso)] > 0:
            rows.append(
                (
                    i,
                    "T2D",
                    iso,
                    "glp1_total",
                    pif * _ow_dict.get(iso, 0.5) / paf_t2d[("T2D", iso)] * 100,
                )
            )
    return rows


def run(job):
    kind, i = job
    return ckd_draw(i) if kind == "CKD" else t2d_draw(i)


if __name__ == "__main__":
    jobs = [("T2D", i) for i in range(A.n_t2d)] + [("CKD", i) for i in range(A.n_ckd)]
    print(f"draws: CKD {A.n_ckd}, T2D {A.n_t2d}; workers {A.workers}", flush=True)
    if A.workers == 1:
        parts = [run(j) for j in jobs]
    else:
        with Pool(A.workers) as p:
            parts = p.map(run, jobs, chunksize=1)
    df = pd.DataFrame(
        [r for part in parts for r in part],
        columns=["draw", "outcome", "iso3", "scenario", "v"],
    )
    df = df.sort_values(["outcome", "draw", "iso3", "scenario"]).reset_index(drop=True)
    df.to_csv(f"results/figure2b_ckd_t2d_draws{A.suffix}.csv", index=False)

    q = (
        df.groupby(["outcome", "iso3", "scenario"])
        .v.quantile([0.5, 0.025, 0.975])
        .unstack()
    )
    q.columns = ["med", "lo", "hi"]
    q = q.reset_index()
    q.to_csv(f"results/figure2b_ckd_t2d_arms_ui{A.suffix}.csv", index=False)

    ratio_rows = []
    for outc in ["CKD", "T2D"]:
        glp1 = df[(df.outcome == outc) & (df.scenario == "glp1_total")]
        for scen in SCENS:
            diet = df[(df.outcome == outc) & (df.scenario == scen)]
            m = glp1[["draw", "iso3", "v"]].merge(
                diet[["draw", "iso3", "v"]], on=["draw", "iso3"], suffixes=("_g", "_d")
            )
            m["ratio"] = m.v_g / m.v_d.replace(0, np.nan)
            m = m.dropna(subset=["ratio"])
            for iso in m.iso3.unique():
                sub = m[m.iso3 == iso]
                if len(sub) < min(10, A.n_t2d):
                    continue
                ratio_rows.append(
                    {
                        "outcome": outc,
                        "iso3": iso,
                        "scenario": scen,
                        "med": sub.ratio.median(),
                        "lo": sub.ratio.quantile(0.025),
                        "hi": sub.ratio.quantile(0.975),
                    }
                )
    ratios = pd.DataFrame(ratio_rows)
    ratios.to_csv(f"results/figure2b_ckd_t2d_ratio_ui{A.suffix}.csv", index=False)
    print(f"arms: {len(q)} rows | ratios: {len(ratios)} rows", flush=True)
    print(ratios[ratios.iso3 == "GLB"].to_string(index=False), flush=True)
