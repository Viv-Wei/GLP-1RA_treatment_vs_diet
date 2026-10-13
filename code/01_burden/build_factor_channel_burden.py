import os

MEASURE = os.environ.get("F1_MEASURE", "Deaths")
SFX = "" if MEASURE == "Deaths" else "_dalys"
import pandas as pd, glob, numpy as np, os

import pandas as _pd_cm

_cm = _pd_cm.read_csv("params/country_map.csv")
ISO = dict(zip(_cm.gbd_location_name, _cm.iso3))
AGE15 = [f"{a}-{a+4} years" for a in range(25, 95, 5)] + ["95+ years"]
RISK2FAC = {
    "Diet low in whole grains": "wholegrains",
    "Diet low in fruits": "fruits",
    "Diet low in vegetables": "vegetables",
    "Diet low in fiber": "fiber",
    "Diet low in fibre": "fiber",
    "Diet low in legumes": "legumes",
    "Diet low in nuts and seeds": "nutsseeds",
    "Diet low in seafood omega-3 fatty acids": "seafood_omega3",
    "Diet low in omega-6 polyunsaturated fatty acids": "omega6",
    "Diet high in red meat": "redmeat",
    "Diet high in processed meat": "procmeat",
    "Diet high in sugar-sweetened beverages": "ssb",
    "Diet high in sodium": "sodium",
    "Diet high in trans fatty acids": "transfat",
}


def load_dir_style():
    DIS = {
        "Ischemic heart disease": "IHD",
        "Ischemic stroke": "IschStroke",
        "Diabetes mellitus type 2": "T2D",
    }
    out = []
    for dis, tag in DIS.items():
        for p in sorted(glob.glob(f"data/gbd/{dis}/*/*/*.csv")):
            risk = p.split(os.sep)[3]
            if risk == "Dietary risks":
                continue
            fac = RISK2FAC.get(risk)
            if fac is None:
                continue
            d = pd.read_csv(p, low_memory=False)
            d["outcome"], d["factor"] = tag, fac
            out.append(d)
    return pd.concat(out, ignore_index=True)


def load_merged_style():
    out = []
    for f, tag in [
        ("data/gbd/queryA_ICH_all_merged.csv", "ICH"),
        ("data/gbd/queryA_SAH_all_merged.csv", "SAH"),
        ("data/gbd/queryA_ckd_3years.csv", "CKD"),
    ]:
        d = pd.read_csv(f, low_memory=False)
        d = d[d.rei_name != "Dietary risks"]
        d["outcome"] = tag
        d["factor"] = d.rei_name.map(RISK2FAC)
        out.append(d[d.factor.notna()])
    return pd.concat(out, ignore_index=True)


def main():
    d = pd.concat([load_dir_style(), load_merged_style()], ignore_index=True)
    d = d[
        (d.year == 2023)
        & d.age_name.isin(AGE15)
        & d.sex_name.isin(["Male", "Female"])
        & d.location_name.isin(ISO)
    ]
    d["iso3"] = d.location_name.map(ISO)

    paf = d[(d.metric_name == "Percent") & (d.measure_name == MEASURE)].rename(
        columns={"val": "paf"}
    )[["iso3", "outcome", "factor", "sex_name", "age_name", "paf"]]
    num = d[(d.metric_name == "Number") & (d.measure_name == MEASURE)].rename(
        columns={"val": "deaths_attr", "upper": "deaths_hi", "lower": "deaths_lo"}
    )[
        [
            "iso3",
            "outcome",
            "factor",
            "sex_name",
            "age_name",
            "deaths_attr",
            "deaths_lo",
            "deaths_hi",
        ]
    ]
    t = paf.merge(num, on=["iso3", "outcome", "factor", "sex_name", "age_name"])

    ch = pd.read_csv("params/glp1_channels_gbd2023.csv")
    print("[channel table columns]", ch.columns.tolist())
    t = t.merge(
        ch[["factor", "outcome", "glp1_channel", "convergence"]],
        on=["factor", "outcome"],
        how="left",
    )
    miss = t[t.glp1_channel.isna()][["factor", "outcome"]].drop_duplicates()
    if len(miss):
        print("[!! no channel label; treated as direct effect only]")
        print(miss.to_string(index=False))
        t["glp1_channel"] = t.glp1_channel.fillna(
            "Direct effect only - structurally unreachable"
        )
        t["convergence"] = t.convergence.fillna(0.0)

    t.to_csv(f"results/factor_burden_by_age{SFX}.csv", index=False)
    fac = (
        t.groupby(["iso3", "outcome", "glp1_channel", "factor"])
        .agg(d=("deaths_attr", "sum"), lo=("deaths_lo", "sum"), hi=("deaths_hi", "sum"))
        .reset_index()
    )
    fac["hw_lo"] = fac.d - fac.lo
    fac["hw_hi"] = fac.hi - fac.d
    import numpy as _np

    g2 = (
        fac.groupby(["iso3", "outcome", "glp1_channel"])
        .apply(
            lambda x: pd.Series(
                dict(
                    deaths=x.d.sum(),
                    deaths_lo=x.d.sum() - _np.sqrt((x.hw_lo**2).sum()),
                    deaths_hi=x.d.sum() + _np.sqrt((x.hw_hi**2).sum()),
                )
            ),
            include_groups=False,
        )
        .reset_index()
    )
    g = (
        t.groupby(["iso3", "outcome", "glp1_channel"])
        .agg(
            deaths_DROP=("deaths_attr", "sum"),
            conv_wt=(
                "deaths_attr",
                lambda s: np.average(
                    t.loc[s.index, "convergence"], weights=s.clip(lower=0) + 1e-9
                ),
            ),
        )
        .reset_index()
    )
    g = g.drop(columns=["deaths_DROP"]).merge(
        g2, on=["iso3", "outcome", "glp1_channel"]
    )
    g[["deaths", "deaths_lo", "deaths_hi"]] = g[
        ["deaths", "deaths_lo", "deaths_hi"]
    ].round(0)
    g.to_csv(f"results/figure1_channel_burden{SFX}.csv", index=False)

    wide = g.pivot_table(
        index=["iso3", "outcome"],
        columns="glp1_channel",
        values="deaths",
        aggfunc="sum",
        fill_value=0,
    )
    chan_cols = list(wide.columns)
    wide["net_total"] = wide[chan_cols].sum(axis=1)
    pos = wide[chan_cols].clip(lower=0).sum(axis=1)
    wide["negative_offset"] = wide[chan_cols].clip(upper=0).sum(axis=1).round(0)
    for c in chan_cols:
        wide[f"{c}%"] = (wide[c].clip(lower=0) / pos * 100).round(1)
    print("\n===== Figure 1 data: attributable deaths (count) and channel shares =====")
    pd.set_option("display.width", 220)
    print(wide.to_string())
    print("\n-> results/figure1_channel_burden.csv")


if __name__ == "__main__":
    main()
