import pandas as pd, numpy as np, argparse


def load_trial(path):
    t = pd.read_csv(path)

    def get(param, trial):
        s = t[(t.param == param) & (t.trial == trial)]
        if not len(s) or pd.isna(s.value.iloc[0]):
            raise SystemExit(f"[missing parameter] {trial}:{param} not filled")
        return float(s.value.iloc[0])

    return get


class MedCurve:
    def __init__(self, sub):
        s = sub.sort_values("exp")
        self.x = s["exp"].to_numpy(float)
        rr = s.rr.to_numpy(float)
        t = float(s.tmrel.iloc[0])
        self.y = rr / np.interp(np.clip(t, self.x.min(), self.x.max()), self.x, rr)

    def rr(self, v):
        return np.maximum(
            np.interp(np.clip(v, self.x.min(), self.x.max()), self.x, self.y), 1.0
        )


def pif_shift(curve, exposure, delta):
    rr0, rr1 = curve.rr(exposure), curve.rr(exposure + delta)
    return float(1.0 - rr1 / rr0)


def main(a):
    get = load_trial(a.trial_params)
    d_sbp = get("sbp_change", "STEP1")
    d_fpg = get("fpg_change", "STEP1")
    retain = 1.0 if not a.adherence_adjust else (1 - 2 / 3 * (1 - a.persist_frac))

    dr = pd.read_csv(a.curves)
    curves = {
        "SBP_IHD": MedCurve(dr[dr.curve_id == "high_sbp_ihd"]),
        "SBP_S": MedCurve(dr[dr.curve_id == "high_sbp_s"]),
        "FPG_IHD": MedCurve(dr[dr.curve_id == "high_glucose_ihd"]),
        "FPG_IS": MedCurve(dr[dr.curve_id == "high_glucose_istroke"]),
    }

    expo = pd.read_csv(a.exposure)
    expo = expo[expo.year_id == 2023]
    import pandas as _cmpd

    _cm = _cmpd.read_csv("params/country_map.csv")
    ISO = dict(zip(_cm.gbd_location_name, _cm.iso3))
    expo["iso3"] = expo.location_name.map(ISO)

    _ow = pd.read_csv("params/overweight_prevalence.csv")
    _ow_dict = dict(zip(_ow.iso3, _ow.overweight_prev))

    rows = []
    for cov in [float(x) for x in a.coverage.split(",")]:
        eff_base = cov * retain
        for (iso, sex, age), g in expo.groupby(["iso3", "sex", "age"]):
            eff = eff_base * _ow_dict.get(iso, 0.5)
            sbp = g[g.rei == "HIGH_SYSTOLIC_BLOOD_PRESSURE"]["mean"]
            fpg = g[g.rei == "HIGH_FASTING_PLASMA_GLUCOSE"]["mean"]
            if not len(sbp) or not len(fpg):
                continue
            sbp, fpg = float(sbp.iloc[0]), float(fpg.iloc[0])
            rows.append(
                dict(
                    coverage=cov,
                    iso3=iso,
                    sex=sex,
                    age=age,
                    sbp0=sbp,
                    fpg0=fpg,
                    pif_ihd_sbp=pif_shift(curves["SBP_IHD"], sbp, eff * d_sbp),
                    pif_str_sbp=pif_shift(curves["SBP_S"], sbp, eff * d_sbp),
                    pif_ihd_fpg=pif_shift(curves["FPG_IHD"], fpg, eff * d_fpg),
                    pif_is_fpg=pif_shift(curves["FPG_IS"], fpg, eff * d_fpg),
                )
            )
    out = pd.DataFrame(rows)
    out.to_csv(a.out, index=False)
    print(f"{len(out):,} rows -> {a.out}")
    print(
        "\n=== Overview: PIF for IHD via SBP at 20% coverage (mean across age and sex, %) ==="
    )
    q = out[out.coverage == 0.20].groupby("iso3").pif_ihd_sbp.mean() * 100
    print(q.round(2).to_string())
    print("\n=== FPG channel, same definition ===")
    q2 = out[out.coverage == 0.20].groupby("iso3").pif_ihd_fpg.mean() * 100
    print(q2.round(2).to_string())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--trial-params", default="params/glp1_trial_params.csv")
    p.add_argument("--curves", default="params/dose_response_gbd2023_gbdexp.csv")
    p.add_argument("--exposure", default="params/exposure_gbd2023.csv")
    p.add_argument("--coverage", default="0.05,0.10,0.20,0.50")
    p.add_argument("--adherence-adjust", action="store_true")
    p.add_argument("--persist-frac", type=float, default=0.5)
    p.add_argument("--out", default="results/glp1_arm_pif.csv")
    main(p.parse_args())
