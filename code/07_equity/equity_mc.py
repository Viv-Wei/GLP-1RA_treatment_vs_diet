import pandas as pd, numpy as np, argparse
from importlib.machinery import SourceFileLoader

dpa = SourceFileLoader(
    "dpa", "code/02_intervention_arms/diet_policy_arm.py"
).load_module()


def main(a):
    rng = np.random.default_rng(a.seed)
    g = pd.read_csv("data/exposure_19country.csv", low_memory=False)
    g = g[g.age >= 22.5]
    nat = (
        g[g.stratum_type == "age_marginal"]
        .groupby(["gdd_var", "iso3"])
        .agg(
            m=("intake_median", "mean"),
            lo=("intake_lo95", "mean"),
            hi=("intake_hi95", "mean"),
        )
    )
    cr = g[g.stratum_type == "cross"]
    st = {}
    for e in [1, 3]:
        st[e] = (
            cr[cr.edu == e]
            .groupby(["gdd_var", "iso3"])
            .agg(
                m=("intake_median", "mean"),
                lo=("intake_lo95", "mean"),
                hi=("intake_hi95", "mean"),
            )
        )
    curves = pd.read_csv("params/dose_response_gbd2023_gbdexp.csv")
    cvt = pd.read_csv("params/cv_by_factor.csv")
    cvmap = dict(zip(cvt.gdd_var, cvt.cv_between * 1.8))
    e23 = pd.read_csv("params/exposure_gbd_male.csv")
    e23 = e23[e23.year == 2023]
    sbpf = "params/sbp_gbd2023_male.csv"

    recs = []
    for d in range(a.draws):
        ratios = {}
        for (gv, iso), r in nat.iterrows():
            sd_n = max((r.hi - r.lo) / (2 * 1.96), 1e-9)
            n_d = max(rng.normal(r.m, sd_n), 1e-9)
            for e in [1, 3]:
                if (gv, iso) in st[e].index:
                    rs = st[e].loc[(gv, iso)]
                    sd_s = max((rs.hi - rs.lo) / (2 * 1.96), 1e-9)
                    s_d = max(rng.normal(rs.m, sd_s), 1e-9)
                    ratios[(gv, iso, e)] = s_d / n_d
        for e in [1, 3]:
            et = e23.copy()
            et["intake_median"] = et.apply(
                lambda r: r.intake_median * ratios.get((r.gdd_var, r.iso3, e), 1.0),
                axis=1,
            )
            paf = dpa.combined_paf(
                et, curves, cvmap, outcomes=("IHD", "IschStroke", "T2D")
            )
            nb = {
                "IHD": dpa.sodium_paf_via_sbp(et, curves, sbpf),
                "IschStroke": dpa.sodium_paf_via_sbp(
                    et, curves, sbpf, curve_id="high_sbp_s"
                ),
                "CKD": dpa.sodium_paf_via_sbp(
                    et, curves, sbpf, curve_id="high_sbp_ckd"
                ),
                "T2D": {iso: 0.0 for iso in set(k[1] for k in paf)},
            }
            for (outc, iso), v in paf.items():
                recs.append(
                    dict(
                        draw=d,
                        edu=e,
                        outcome=outc,
                        iso3=iso,
                        paf=1 - (1 - v) * (1 - nb[outc].get(iso, 0.0)),
                    )
                )
            for iso, ckd_paf in nb["CKD"].items():
                recs.append(dict(draw=d, edu=e, outcome="CKD", iso3=iso, paf=ckd_paf))
        if (d + 1) % 10 == 0:
            print(f"{d+1}/{a.draws}")
    x = pd.DataFrame(recs)
    w = x.pivot_table(index=["draw", "outcome", "iso3"], columns="edu", values="paf")
    w["ratio"] = w[1] / w[3]
    q = (
        w.groupby(["outcome", "iso3"])
        .ratio.quantile([0.5, 0.025, 0.975])
        .unstack()
        .round(3)
    )
    q.columns = ["med", "lo", "hi"]
    q.to_csv("results/equity_ratio_ui.csv")
    print(q.to_string())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--draws", type=int, default=200)
    p.add_argument("--seed", type=int, default=42)
    main(p.parse_args())
