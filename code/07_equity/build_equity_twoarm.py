import argparse, os, csv
import numpy as np, pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--root", default=".")
a = p.parse_args()
R = lambda f: os.path.join(a.root, f)

ISO_W = {
    "AUS": "Australia",
    "BGD": "Bangladesh",
    "BGR": "Bulgaria",
    "BRA": "Brazil",
    "CAN": "Canada",
    "CHN": "China",
    "DEU": "Germany",
    "EGY": "Egypt",
    "GBR": "United Kingdom of Great Britain and Northern Ireland",
    "IND": "India",
    "ISR": "Israel",
    "JPN": "Japan",
    "MEX": "Mexico",
    "NGA": "Nigeria",
    "NZL": "New Zealand",
    "ROU": "Romania",
    "TUR": "Turkey",
    "USA": "United States of America",
    "ZAF": "South Africa",
}
ALT = {
    "GBR": ["United Kingdom"],
    "TUR": ["Türkiye", "Turkiye"],
    "USA": ["United States"],
    "EGY": ["Egypt, Arab Rep."],
}
EDU3 = {
    "No Education": 1,
    "Incomplete Primary": 1,
    "Primary": 1,
    "Lower Secondary": 2,
    "Upper Secondary": 2,
    "Post Secondary": 3,
    "Short Post Secondary": 3,
    "Bachelor": 3,
    "Master and higher": 3,
}
LAB = {1: "Low education", 2: "Medium education", 3: "High education"}
PRICES = {
    "Branded": "afford_brand_pct",
    "Generic $140/yr": "afford_generic_hi_pct",
    "Generic $28/yr": "afford_generic_lo_pct",
}
OUT_W = {
    "IHD": "Ischaemic heart disease",
    "IschStroke": "Ischaemic stroke",
    "CKD": "Chronic kidney disease",
}
THRESH = 10.0

w = pd.read_csv(R("data/wittgenstein/wcde_clean.csv"))
yr = int(w.Year.unique()[np.argmin(np.abs(w.Year.unique() - 2023))])
w = w[(w.Year == yr) & w.Education.isin(EDU3)].copy()
w["a0"] = pd.to_numeric(w.Age.astype(str).str.extract(r"^(\d+)")[0], errors="coerce")
w = w[w.a0 >= 25]
share = {}
for iso, nm in ISO_W.items():
    x = w[w.Area.isin([nm] + ALT.get(iso, []))]
    assert len(x), f"Wittgenstein: no rows for {iso} ({nm})"
    g = x.assign(e=x.Education.map(EDU3)).groupby("e").Population.sum()
    share[iso] = (g / g.sum()).reindex([1, 2, 3]).fillna(0)
S = pd.DataFrame(share).T

paf = pd.read_csv(R("results/equity_paf_by_edu.csv"))
paf = paf[paf.outcome.isin(OUT_W)]
j = pd.read_csv(R("results/supp_table2_joint.csv"))
dcol = next(c for c in j.columns if "death" in c.lower() and "rate" not in c.lower())
ocol = next(c for c in j.columns if c.lower() in ("outcome", "cause"))
icol = next(c for c in j.columns if c.lower() in ("iso3", "iso"))
wts = {}
for iso in ISO_W:
    jj = j[j[icol] == iso].set_index(ocol)[dcol]
    v = pd.Series({o: float(jj.get(o, jj.get(OUT_W[o], np.nan))) for o in OUT_W}).clip(
        lower=0
    )
    assert v.notna().all(), f"joint burden missing for {iso}: {v.to_dict()}"
    wts[iso] = v / v.sum()


def diet_alloc(iso):
    s = S.loc[iso]
    out = pd.Series(0.0, index=[1, 2, 3])
    for o, wo in wts[iso].items():
        q = (
            paf[(paf.iso3 == iso) & (paf.outcome == o)]
            .set_index("edu")
            .paf.reindex([1, 2, 3])
        )
        if q.isna().any() or (s * q).sum() <= 0:
            out += wo * s
        else:
            out += wo * (s * q) / (s * q).sum()
    return out


def paf_ratio(iso):
    s = S.loc[iso]
    r = pd.Series(0.0, index=[1, 2, 3])
    for o, wo in wts[iso].items():
        q = (
            paf[(paf.iso3 == iso) & (paf.outcome == o)]
            .set_index("edu")
            .paf.reindex([1, 2, 3])
        )
        r += wo * (q / (s * q).sum() if q.notna().all() and (s * q).sum() > 0 else 1.0)
    return r


aff = pd.read_csv(R("results/affordability_matrix.csv"))


def afford_frac(iso, col, lam=1.0):
    ok = (
        aff[aff.iso3 == iso]
        .set_index("quintile")[col]
        .reindex(["Q1", "Q2", "Q3", "Q4", "Q5"])
        <= THRESH
    )
    s = S.loc[iso]
    cum = np.concatenate([[0], np.cumsum(s.values)])
    out = {}
    for k, e in enumerate([1, 2, 3]):
        lo_, hi_ = cum[k], cum[k + 1]
        if hi_ - lo_ <= 0:
            out[e] = 0.0
            continue
        cover = sum(
            max(0, min(hi_, (q + 1) * 0.2) - max(lo_, q * 0.2)) * ok.iloc[q]
            for q in range(5)
        )
        out[e] = cover / (hi_ - lo_)
    rank = pd.Series(out)
    avg = ok.mean()
    return lam * rank + (1 - lam) * avg


def mean_afford_pct(iso, col):
    v = (
        aff[aff.iso3 == iso]
        .set_index("quintile")[col]
        .reindex(["Q1", "Q2", "Q3", "Q4", "Q5"])
        .values
    )
    s = S.loc[iso]
    cum = np.concatenate([[0], np.cumsum(s.values)])
    out = {}
    for k, e in enumerate([1, 2, 3]):
        wq = [
            max(0, min(cum[k + 1], (q + 1) * 0.2) - max(cum[k], q * 0.2))
            for q in range(5)
        ]
        out[e] = float(np.dot(wq, v) / sum(wq)) if sum(wq) > 0 else np.nan
    return out


t9 = pd.read_csv(R("results/itt_1012/tableS9_deaths_itt.csv"))
t9 = t9[t9.coverage == "20% coverage"].pivot_table(
    index="iso3", columns="arm", values="med"
)
DIET, GLP = t9["Dietary policy (feasible)"], t9["GLP-1RA treatment"]

rows = []
for price, col in PRICES.items():
    for iso in sorted(ISO_W):
        da, pr, af, ap = (
            diet_alloc(iso),
            paf_ratio(iso),
            afford_frac(iso, col),
            mean_afford_pct(iso, col),
        )
        for e in [1, 2, 3]:
            d_e = DIET[iso] * da[e] / 1e3
            g_e = GLP[iso] * S.loc[iso, e] * af[e] / 1e3
            rows.append(
                [
                    price,
                    iso,
                    LAB[e],
                    round(S.loc[iso, e], 3),
                    round(pr[e], 2),
                    round(ap[e], 1),
                    round(af[e], 3),
                    round(d_e, 1),
                    round(g_e, 1),
                    round(d_e / g_e, 2) if g_e > 0 else None,
                ]
            )
cols = [
    "Price scenario",
    "ISO3",
    "Education group",
    "Population share",
    "PAF ratio",
    "Affordability (% income)",
    "Affordable fraction",
    "Dietary policy averted (thousands)",
    "GLP-1RA averted (thousands)",
    "Diet-to-treatment ratio",
]
df = pd.DataFrame(rows, columns=cols)
os.makedirs(R("submission/tables"), exist_ok=True)
with open(R("tables/tableS_equity_twoarm.csv"), "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(
        [
            "Supplementary Table. Deaths averted over 30 years by dietary policy and GLP-1RA treatment, by education group and drug price."
        ]
    )
    wr.writerow(
        [
            f"Education structure: Wittgenstein Centre {yr}, adults aged >=25 years. Dietary policy deaths averted allocated "
            "by population share x education-specific PAF (IHD, ischaemic stroke, CKD, weighted by diet-attributable deaths). "
            "Education groups placed on the population percentile scale in rank order; income quintiles affordable if "
            "annual cost <=10% of per-capita income. Treatment deaths averted = national x population share x affordable fraction."
        ]
    )
    wr.writerow([])
    wr.writerow(cols)
    for r in df.itertuples(index=False):
        wr.writerow(list(r))

sens = []
for lam, lab in [
    (1.0, "Education and income perfectly ranked (main analysis)"),
    (0.5, "Partial link (50%)"),
    (0.0, "Education unrelated to income"),
]:
    dl = gl = 0.0
    none = []
    for iso in sorted(ISO_W):
        af = afford_frac(iso, PRICES["Generic $140/yr"], lam)
        dl += DIET[iso] * diet_alloc(iso)[1]
        g = GLP[iso] * S.loc[iso, 1] * af[1]
        gl += g
        if g <= 0:
            none.append(iso)
    sens.append(
        {
            "Education-income link": lab,
            "Dietary policy deaths averted (thousands)": round(dl / 1e3),
            "GLP-1RA treatment deaths averted (thousands)": round(gl / 1e3),
            "Diet-to-treatment ratio": round(dl / gl, 2),
            "Countries with no treatment benefit": f"{len(none)}/19"
            + (f" ({', '.join(none)})" if none else ""),
        }
    )
pd.DataFrame(sens).to_csv(
    R("tables/tableS_equity_mapping_sensitivity.csv"), index=False
)

print(f"Wittgenstein year {yr}; education shares (low/medium/high), first rows:")
print(S.round(2).head(6).to_string())
pool = df.groupby(["Price scenario", "Education group"])[
    ["Dietary policy averted (thousands)", "GLP-1RA averted (thousands)"]
].sum()
print("\nPooled over 19 countries (thousands):\n" + pool.round(0).to_string())
print(
    "\nEducation-income link (low education, generic US$140):\n"
    + pd.DataFrame(sens).to_string(index=False)
)
chk = (
    df[df["Price scenario"] == "Branded"]
    .groupby("ISO3")["Dietary policy averted (thousands)"]
    .sum()
    - DIET / 1e3
)
print(
    f"\nCheck: dietary policy allocation sums to national total (max abs difference, thousands): {chk.abs().max():.3f}"
)
