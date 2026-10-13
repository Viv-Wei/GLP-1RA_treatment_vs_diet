import argparse, os
import numpy as np, pandas as pd, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--root", default=".")
p.add_argument("--out", default="figures/figure5_robustness")
a = p.parse_args()
R = lambda f: os.path.join(a.root, f)
lead = lambda df, col="ratio": int((df.loc[df.iso3 != "GLB", col] > 1).sum())

t10 = pd.read_csv(R("results/itt_1012/tableS10_ranking_threshold_itt.csv")).set_index(
    "iso3"
)
sd = pd.read_csv(R("results/itt_1012/sens_deterministic_itt.csv"))
b30 = pd.read_csv(R("results/itt_1012/sens_bmi30_itt.csv"))
dc = pd.read_csv(R("results/itt_1012/sens_discount_itt.csv"))
sw = pd.read_csv(R("results/sensitivity_sweep_19c.csv"))
sw["ratio"] = sw.diet_30yr_k / sw.glp1_30yr_k
st = pd.read_csv(
    R("tables/tableS_lifetable_structural_sensitivity.csv"), skiprows=3
).set_index("Scenario")
_fsc = R("results/itt_1012/sens_structural_countries_itt.csv")
stc = pd.read_csv(_fsc) if os.path.exists(_fsc) else None


def struct(scn):
    if stc is None:
        return st.loc[scn].iloc[3], None
    x = stc[stc.scenario == scn].assign(ratio=lambda t: t.diet / t.glp1)
    return float(x.loc[x.iso3 == "GLB", "ratio"].iloc[0]), lead(x)


po = pd.read_csv(R("tables/tableS_population_sensitivity_v2.csv"), skiprows=3).rename(
    columns={"Country": "iso3"}
)


def det(scn, outc="all"):
    x = sd[(sd.scenario == scn) & (sd.outcomes == outc)]
    return float(x.loc[x.iso3 == "GLB", "ratio"].iloc[0]), lead(x)


def sweep(pm, rt, rd):
    x = sw[
        np.isclose(sw.policy_maintain, pm)
        & np.isclose(sw.glp1_retain, rt)
        & (sw.ramp_diet == rd)
    ]
    return float(x.loc[x.iso3 == "GLB", "ratio"].iloc[0]), lead(x)


def disc(rate):
    x = dc[np.isclose(dc.rate, rate)]
    return float(x.loc[x.iso3 == "GLB", "ratio"].iloc[0]), lead(x)


base = sd[(sd.scenario == "base_0.20") & (sd.outcomes == "all")].copy()
base["r_itt_off"] = base.diet / (base.glp1 * 0.667 / 0.76)
g = t10.loc["GLB"]
main = (
    g.threshold_med / 20,
    g.threshold_lo / 20,
    g.threshold_hi / 20,
    lead(t10.reset_index().assign(ratio=t10.diet_med.values / t10.glp1_med.values)),
)

ROWS = [
    ("Model structure", "Incidence share of treatment effect 0%", *det("split_0")),
    ("Model structure", "Incidence share of treatment effect 100%", *det("split_1.0")),
    (
        "Model structure",
        "Same-cycle case fatality for 50% of new cases",
        *struct("Same-cycle CFR for 50% of new cases"),
    ),
    (
        "Model structure",
        "Background mortality in diseased ×1.5",
        *struct("Background mortality 1.5x for diseased"),
    ),
    (
        "Model structure",
        "Background mortality in diseased ×2.0",
        *struct("Background mortality 2.0x for diseased"),
    ),
    (
        "Model structure",
        "Background mortality in diseased ×3.0",
        *struct("Background mortality 3.0x for diseased"),
    ),
    (
        "Treatment effect and eligibility",
        "Cardiovascular HR 0.74 (SUSTAIN-6)",
        *det("cvhr_0.74"),
    ),
    (
        "Treatment effect and eligibility",
        "Cardiovascular HR 0.86 (SOUL)",
        *det("cvhr_0.86"),
    ),
    (
        "Treatment effect and eligibility",
        "CKD HR 0.66 (FLOW lower 95% CI)",
        *det("ckdhr_0.66"),
    ),
    (
        "Treatment effect and eligibility",
        "CKD HR 0.88 (FLOW upper 95% CI)",
        *det("ckdhr_0.88"),
    ),
    (
        "Treatment effect and eligibility",
        "CKD excluded (IHD and IS only)",
        *det("base_0.20", "noCKD"),
    ),
    (
        "Treatment effect and eligibility",
        "Eligibility restricted to BMI ≥30",
        float(b30.loc[b30.iso3 == "GLB", "ratio"].iloc[0]),
        lead(b30),
    ),
    ("Treatment persistence", "Retention 0.37 (real-world 0.33)", *sweep(1.0, 0.37, 8)),
    ("Treatment persistence", "Retention 0.57 (real-world 0.50)", *sweep(1.0, 0.57, 8)),
    (
        "Treatment persistence",
        "Retention 1.13 (full persistence)",
        *sweep(1.0, 1.13, 8),
    ),
    (
        "Treatment persistence",
        "Retention 0.667 (no ITT correction)",
        float(base.loc[base.iso3 == "GLB", "r_itt_off"].iloc[0]),
        lead(base, "r_itt_off"),
    ),
    (
        "Dietary policy implementation",
        "Policy effect decays to 80%",
        *sweep(0.8, 0.76, 8),
    ),
    (
        "Dietary policy implementation",
        "Policy effect decays to 60%",
        *sweep(0.6, 0.76, 8),
    ),
    (
        "Dietary policy implementation",
        "Policy effect decays to 40%",
        *sweep(0.4, 0.76, 8),
    ),
    ("Dietary policy implementation", "Ramp-up over 12 years", *sweep(1.0, 0.76, 12)),
    ("Dietary policy implementation", "Ramp-up over 16 years", *sweep(1.0, 0.76, 16)),
    (
        "Dietary policy implementation",
        "Decay 40%, 16-year ramp-up, retention 1.13",
        *sweep(0.4, 1.13, 16),
    ),
    (
        "Population and time",
        "2053 population structure",
        float(po.loc[po.iso3 == "GLB", "Ratio 2053"].iloc[0]),
        int((po.loc[po.iso3 != "GLB", "Ratio 2053"] > 1).sum()),
    ),
    ("Population and time", "Discounting 3% a year", *disc(0.03)),
    ("Population and time", "Discounting 5% a year", *disc(0.05)),
]


def _newgroup(g, label):
    if label.startswith("Eligibility restricted"):
        return "Eligible population"
    if "no ITT correction" in label:
        return "Correction for discontinuation in trials"
    if g == "Treatment effect and eligibility":
        return "Trial effect estimates"
    return g


SECTIONS = {
    "Conditions": [
        "Treatment persistence",
        "Dietary policy implementation",
        "Eligible population",
    ],
    "Sensitivity analyses": [
        "Model structure",
        "Trial effect estimates",
        "Correction for discontinuation in trials",
        "Population and time",
    ],
}
GORDER = ["Main analysis"] + [g for gs in SECTIONS.values() for g in gs]
SECTION_OF = {g: s for s, gs in SECTIONS.items() for g in gs}
ROWS = [(_newgroup(r[0], r[1]),) + tuple(r[1:]) for r in ROWS]
ROWS = sorted(ROWS, key=lambda r: GORDER.index(r[0]))
tab = pd.DataFrame(ROWS, columns=["group", "label", "ratio", "n_diet_leads"])
tab.loc[-1] = [
    "Main analysis",
    "Main analysis (Monte Carlo median, 95% UI)",
    main[0],
    main[3],
]
tab = tab.sort_index().reset_index(drop=True)
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
tab.assign(threshold_pct=20 * tab.ratio).round(3).to_csv(
    a.out + "_data.csv", index=False
)

GCOL = {
    "Main analysis": "#1a1a1a",
    "Treatment persistence": "#B5651D",
    "Dietary policy implementation": "#00A087",
    "Eligible population": "#E64B35",
    "Model structure": "#7F7F7F",
    "Trial effect estimates": "#8491B4",
    "Correction for discontinuation in trials": "#A67C52",
    "Population and time": "#3C5488",
}
items, y, cur_sec = [], 0, None
for grp, sub in tab.groupby("group", sort=False):
    if grp != "Main analysis":
        if SECTION_OF[grp] != cur_sec:
            cur_sec = SECTION_OF[grp]
            y += 0.9
            items.append(("section", cur_sec, y))
            y += 0.6
        y += 0.5
        items.append(("header", grp, y))
        y += 1
    for r in sub.itertuples():
        items.append(("row", r, y))
        y += 1
H = y

plt.rcParams.update(
    {
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "pdf.fonttype": 42,
    }
)
fig = plt.figure(figsize=(7.4, 7.0))
ax = fig.add_axes([0.40, 0.07, 0.40, 0.84])
XL = (0.3, 8)
ax.set_xscale("log")
ax.set_xlim(*XL)
ax.set_ylim(H - 0.3, -1.0)
ax.axvspan(XL[0], 1, color="#E64B35", alpha=0.06, lw=0)
ax.axvspan(main[1], main[2], color="grey", alpha=0.12, lw=0)
ax.axvline(1, color="black", lw=0.8)
ax.axvline(main[0], color="grey", lw=0.7, ls=(0, (3, 2)))
ax.text(
    0.97,
    H - 0.55,
    "GLP-1RA treatment\naverts more",
    ha="right",
    va="bottom",
    fontsize=6.3,
    color="#b23a29",
    style="italic",
)
ax.text(
    7.7,
    H - 0.55,
    "Dietary policy\naverts more",
    ha="right",
    va="bottom",
    fontsize=6.3,
    color="#00795f",
    style="italic",
)

for kind, obj, yy in items:
    if kind == "section":
        ty = ax.transData.transform((1, yy))[1] / fig.dpi / fig.get_size_inches()[1]
        fig.text(
            0.015, ty, obj, fontsize=8, fontweight="bold", color="black", va="center"
        )
        fig.add_artist(
            plt.Line2D(
                [0.015, 0.985],
                [ty - 0.012] * 2,
                color="black",
                lw=0.6,
                transform=fig.transFigure,
            )
        )
        continue
    if kind == "header":
        fig.text(
            0.015,
            ax.transData.transform((1, yy))[1] / fig.dpi / fig.get_size_inches()[1],
            obj,
            fontsize=7.2,
            fontweight="bold",
            color=GCOL[obj],
            va="center",
        )
        continue
    r = obj
    col = GCOL[r.group]
    is_main = r.group == "Main analysis"
    if is_main:
        ax.plot([main[1], main[2]], [yy, yy], color=col, lw=1.2, zorder=3)
        ax.scatter(r.ratio, yy, marker="D", s=40, color=col, zorder=5)
    else:
        ax.plot([main[0], r.ratio], [yy, yy], color=col, lw=0.8, alpha=0.45, zorder=2)
        ax.scatter(r.ratio, yy, s=26, color=col, edgecolor="white", lw=0.5, zorder=4)
    ty = ax.transData.transform((1, yy))[1] / fig.dpi / fig.get_size_inches()[1]
    fig.text(
        0.03,
        ty,
        r.label,
        fontsize=6.6,
        va="center",
        fontweight="bold" if is_main else "normal",
    )
    txt = f"{r.ratio:.2f}" + (f" ({main[1]:.2f}–{main[2]:.2f})" if is_main else "")
    fig.text(
        0.812,
        ty,
        txt,
        fontsize=6.6,
        va="center",
        fontweight="bold" if is_main else "normal",
    )
    fig.text(
        0.955,
        ty,
        "–" if pd.isna(r.n_diet_leads) else f"{int(r.n_diet_leads)}",
        fontsize=6.6,
        va="center",
        ha="center",
        fontweight="bold" if is_main else "normal",
    )

top = ax.transAxes.transform((0, 1))[1] / fig.dpi / fig.get_size_inches()[1]
fig.text(
    0.812, top + 0.008, "Global ratio", fontsize=6.6, fontweight="bold", va="bottom"
)
fig.text(
    0.955,
    top + 0.008,
    "Diet leads\n(of 19)",
    fontsize=6.6,
    fontweight="bold",
    va="bottom",
    ha="center",
)
fig.text(0.015, top + 0.008, "Analysis", fontsize=7, fontweight="bold", va="bottom")

ticks = [0.3, 0.5, 1, 2, 3, 5]
ax.set_xticks(ticks)
ax.set_xticklabels([f"{t:g}" for t in ticks])
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set_xlabel(
    "Global ratio of deaths averted, dietary policy / GLP-1RA treatment\nat 20% coverage (log scale)"
)
ax.set_yticks([])
ax.spines[["left", "right"]].set_visible(False)
ax2 = ax.secondary_xaxis("top", functions=(lambda r: 20 * r, lambda t: t / 20))
tt = [10, 20, 50, 100]
ax2.set_xticks(tt)
ax2.set_xticklabels([f"{t}%" for t in tt])
ax2.xaxis.set_minor_formatter(plt.NullFormatter())
ax2.set_xlabel(
    "Global coverage needed to match dietary policy (% of adults with BMI ≥25)",
    labelpad=4,
)
fig.savefig(a.out + ".pdf")
fig.savefig(a.out + ".png", dpi=400)
print(tab.round(2).to_string(index=False))
print("saved:", a.out + ".pdf / .png / _data.csv")
