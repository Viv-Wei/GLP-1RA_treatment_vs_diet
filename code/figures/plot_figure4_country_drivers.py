import argparse, os
import numpy as np, pandas as pd, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

p = argparse.ArgumentParser()
p.add_argument("--root", default=".")
p.add_argument("--xlsx", default="tables/Supplementary_Tables_final1012_OW.xlsx")
p.add_argument("--out", default="figures/figure4_country_drivers")
a = p.parse_args()
R = lambda f: os.path.join(a.root, f)
tab = lambda n: pd.read_excel(R(a.xlsx), f"Supplementary Table {n}", header=1).iloc[:-1]

t6, t10, t11 = tab(6), tab(10), tab(11)
name2iso = dict(zip(t6.Country, t6.ISO3))
t10["iso3"] = t10.Country.map(name2iso)
thr = t10.set_index("iso3")["Coverage needed to match dietary policy, median (%)"]
ORDER = ["GLB"] + thr.drop("GLB").sort_values().index.tolist()

PKG = ["Sodium (high)", "Whole grains (low)", "SSB (high)"]
cv = t6.assign(v=t6["Diet-attributable CV deaths"].clip(lower=0))
pkg = (
    cv[cv["Dietary factor"].isin(PKG)].groupby("ISO3").v.sum()
    / cv.groupby("ISO3").v.sum()
    * 100
)

fo = R("params/overweight_prevalence.csv")
if os.path.exists(fo):
    ow = pd.read_csv(fo).set_index("iso3").overweight_prev * 100
else:
    print(
        "WARNING: params/overweight_prevalence.csv not found; using Table 22 (GBD 2023, adults >=25)"
    )
    t22 = tab(22)
    ow = pd.Series(t22["BMI≥25 prev"].values * 100, index=t22.Country.map(name2iso))

d = pd.DataFrame({"thr": thr, "pkg": pkg, "ow": ow}).reindex(ORDER)
d["idx"] = d.pkg / d.ow
cty = d.drop("GLB", errors="ignore").dropna()
b1, b0 = np.polyfit(np.log10(cty.idx), np.log10(cty.thr), 1)
r = np.corrcoef(np.log10(cty.idx), np.log10(cty.thr))[0, 1]

FAC = [
    "Sodium",
    "Whole grains",
    "SSB",
    "Fruits",
    "Vegetables",
    "Legumes",
    "Nuts & seeds",
    "Fibre",
    "Omega-6",
    "Seafood omega-3",
]
sh = t11.pivot_table(
    index="Dietary factor", columns="ISO3", values="Share of country CVD benefit (%)"
).reindex(FAC)[ORDER]
os.makedirs(os.path.dirname(R(a.out)) or ".", exist_ok=True)
d.round(3).to_csv(R(a.out) + "_data.csv")
sh.round(1).to_csv(R(a.out) + "_priorities.csv")

plt.rcParams.update(
    {
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.2,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "pdf.fonttype": 42,
    }
)
BINS = [
    (0, 20, "≤20%", "#E64B35"),
    (20, 50, "20–50%", "#F39B7F"),
    (50, 100, "50–100%", "#4DBBD5"),
    (100, 1e9, ">100%", "#3C5488"),
]
tcol = lambda v: next(
    c for lo, hi, _, c in BINS if lo < v <= hi or (lo == 0 and v <= hi)
)
fig = plt.figure(figsize=(7.4, 6.4))
axA = fig.add_axes([0.075, 0.58, 0.38, 0.36])
axB = fig.add_axes([0.585, 0.58, 0.38, 0.36])
axC = fig.add_axes([0.135, 0.10, 0.83, 0.36])

OFF_A = {
    "EGY": (5, -6),
    "GBR": (-17, 3),
    "AUS": (4, -8),
    "NZL": (-17, 2),
    "DEU": (4, -7),
    "CAN": (4, -8),
    "USA": (4, 2),
    "BRA": (4, -2),
    "ROU": (-20, -7),
    "BGR": (4, -2),
}
OFF_B = {
    "AUS": (-19, 4),
    "NZL": (-20, -2),
    "GBR": (-20, -8),
    "DEU": (2, -9),
    "USA": (4, 1),
    "CAN": (4, -6),
    "TUR": (3, -8),
    "ISR": (4, 0),
    "BRA": (3, 3),
    "ROU": (-20, 2),
}
for iso, rr in d.iterrows():
    glb = iso == "GLB"
    axA.scatter(
        rr.pkg,
        rr.ow,
        s=60 if glb else 34,
        marker="*" if glb else "o",
        color="black" if glb else tcol(rr.thr),
        edgecolor="white",
        lw=0.5,
        zorder=3,
    )
    axA.annotate(
        "Global" if glb else iso,
        (rr.pkg, rr.ow),
        xytext=OFF_A.get(iso, (3, 2)),
        textcoords="offset points",
        fontsize=6,
        fontweight="bold" if glb else "normal",
    )
axA.text(
    1,
    91,
    "Package targets little, many eligible:\nlow coverage needed",
    fontsize=6.2,
    color="grey",
    va="top",
)
axA.text(
    54,
    16,
    "Package targets much, few eligible:\nhigh coverage needed",
    fontsize=6.2,
    color="grey",
    ha="right",
    va="bottom",
)
axA.set_xlim(0, 55)
axA.set_ylim(15, 92)
axA.set_xlabel(
    "Share of diet-attributable cardiovascular deaths\ntargeted by the feasible package (%)"
)
axA.set_ylabel("Overweight prevalence (BMI ≥25, %)")
axA.spines[["top", "right"]].set_visible(False)
axA.legend(
    handles=[
        Line2D([], [], marker="o", ls="", color=c, ms=5, label=l)
        for lo_, hi_, l, c in BINS
        if ((cty.thr > lo_) & (cty.thr <= hi_)).any()
        or (lo_ == 0 and (cty.thr <= hi_).any())
    ],
    title="Coverage needed to\nmatch dietary policy",
    title_fontsize=6.3,
    fontsize=6.2,
    loc="lower left",
    frameon=False,
    handletextpad=0.2,
    borderaxespad=0.2,
)

xs = np.logspace(np.log10(cty.idx.min() * 0.8), np.log10(cty.idx.max() * 1.25), 50)
axB.plot(
    xs, 10 ** (b0 + b1 * np.log10(xs)), color="grey", lw=0.8, ls=(0, (3, 2)), zorder=1
)
for iso, rr in d.iterrows():
    glb = iso == "GLB"
    axB.scatter(
        rr.idx,
        rr.thr,
        s=60 if glb else 34,
        marker="*" if glb else "o",
        color="black" if glb else tcol(rr.thr),
        edgecolor="white",
        lw=0.5,
        zorder=3,
    )
    axB.annotate(
        "Global" if glb else iso,
        (rr.idx, rr.thr),
        xytext=OFF_B.get(iso, (3, 2)),
        textcoords="offset points",
        fontsize=6,
        fontweight="bold" if glb else "normal",
    )
for v, lab in [(20, "Near-term (20%)"), (100, "Full coverage")]:
    axB.axhline(v, color="grey", lw=0.5, ls=":")
    axB.text(xs.min(), v * 1.05, lab, fontsize=6, color="grey")
axB.text(
    0.97,
    0.05,
    f"r = {r:.2f} (log scale, 19 countries)",
    transform=axB.transAxes,
    ha="right",
    fontsize=6.5,
)
axB.set_xscale("log")
axB.set_yscale("log")
axB.set_xticks([0.1, 0.2, 0.5, 1, 2])
axB.set_xticklabels(["0.1", "0.2", "0.5", "1", "2"])
axB.set_yticks([1, 2, 5, 10, 20, 50, 100])
axB.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"])
axB.minorticks_off()
axB.set_xlabel("Targeted share / overweight prevalence")
axB.set_ylabel("Coverage needed to match dietary policy\n(% of adults with BMI ≥25)")
axB.spines[["top", "right"]].set_visible(False)

X = np.arange(len(ORDER))
Y = np.arange(len(FAC))
axC.axhspan(-0.5, 2.5, color="#F3D27A", alpha=0.18, lw=0)
SC = 7.0
for j, iso in enumerate(ORDER):
    col_vals = sh[iso]
    top = col_vals.idxmax()
    for i, f in enumerate(FAC):
        v = col_vals[f]
        if pd.isna(v) or v <= 0:
            continue
        c = "#B5651D" if i < 3 else "#00A087"
        axC.scatter(
            j,
            i,
            s=v * SC,
            color=c,
            alpha=0.85,
            edgecolor="black" if f == top else "white",
            lw=1.0 if f == top else 0.4,
            zorder=3,
        )
        if v >= 15:
            axC.text(
                j,
                i,
                f"{v:.0f}",
                ha="center",
                va="center",
                fontsize=5.4,
                color="white",
                zorder=4,
            )
axC.axvline(0.5, color="grey", lw=0.5)
axC.set_xlim(-0.6, len(ORDER) - 0.4)
axC.set_ylim(len(FAC) - 0.5, -0.5)
axC.set_xticks(X)
axC.set_xticklabels(["Global" if i == "GLB" else i for i in ORDER], fontsize=6.6)
axC.set_yticks(Y)
axC.set_yticklabels(FAC, fontsize=6.8)
for t in axC.get_yticklabels()[:3]:
    t.set_color("#8a4a12")
    t.set_fontweight("bold")
axC.tick_params(length=0)
axC.grid(color="#eeeeee", lw=0.5, zorder=0)
axC.set_axisbelow(True)
[s.set_visible(False) for s in axC.spines.values()]
leg = [
    Line2D(
        [],
        [],
        marker="o",
        ls="",
        color="grey",
        ms=np.sqrt(v * SC),
        alpha=0.7,
        label=f"{v}%",
    )
    for v in [10, 30, 50]
]
leg += [
    Line2D(
        [],
        [],
        marker="o",
        ls="",
        mfc="#B5651D",
        mec="white",
        ms=6,
        label="In feasible package",
    ),
    Line2D(
        [],
        [],
        marker="o",
        ls="",
        mfc="#00A087",
        mec="white",
        ms=6,
        label="Not in package",
    ),
    Line2D(
        [],
        [],
        marker="o",
        ls="",
        mfc="white",
        mec="black",
        mew=1.0,
        ms=6,
        label="Top priority",
    ),
]
axC.legend(
    handles=leg,
    loc="upper center",
    bbox_to_anchor=(0.45, -0.065),
    ncol=6,
    fontsize=6.3,
    frameon=False,
    title="Share of the country's cardiovascular benefit from shifting each factor 20% toward its TMREL",
    title_fontsize=6.4,
    handletextpad=0.6,
    columnspacing=2.0,
    labelspacing=1.2,
)
for ax, tag, dx in [(axA, "a", -0.17), (axB, "b", -0.17), (axC, "c", -0.12)]:
    ax.text(dx, 1.05, tag, transform=ax.transAxes, fontsize=12, fontweight="bold")
fig.savefig(R(a.out) + ".pdf")
fig.savefig(R(a.out) + ".png", dpi=400)
print(f"saved: {R(a.out)} | r = {r:.2f}, slope = {b1:.2f}")
