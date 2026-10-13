import argparse, os
import numpy as np, pandas as pd, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

p = argparse.ArgumentParser()
p.add_argument("--root", default=".")
p.add_argument("--xlsx", default="params/Supplementary_Tables_final1013_HR.xlsx")
p.add_argument("--out", default="figures/suppfig4_static")
a = p.parse_args()
R = lambda f: os.path.join(a.root, f)
tab = lambda n: pd.read_excel(R(a.xlsx), f"Supplementary Table {n}", header=1).iloc[:-1]

t8, t10 = tab(8), tab(10)
for c in t8.columns[4:]:
    t8[c] = pd.to_numeric(t8[c], errors="coerce")
name2iso = dict(zip(t8.Country, t8.ISO3))
t10["iso3"] = t10.Country.map(name2iso)
ORDER = ["GLB"] + t10[t10.iso3 != "GLB"].sort_values(t10.columns[6]).iso3.tolist()
OUTC = [
    ("IHD", "Ischaemic heart disease"),
    ("Ischaemic stroke", "Ischaemic stroke"),
    ("CKD", "Chronic kidney disease"),
    ("T2D", "Type 2 diabetes (incidence)"),
]
SC = [
    ("Sodium reduction 30%", "Sodium reduction 30%", "o", "#3C5488", 0.22),
    ("Feasible policy package", "Feasible package", "D", "#00A087", 0.0),
    ("All factors 20% toward TMREL", "Aspirational scenario", "^", "#7CC9B5", -0.22),
]
GRP = {
    "GLB": "#1a1a1a",
    **{k: "#00A087" for k in "BGR ROU IND MEX ZAF CHN TUR BRA NZL CAN AUS ISR".split()},
    **{k: "#3C5488" for k in "EGY BGD NGA".split()},
    **{k: "#E64B35" for k in "USA DEU GBR JPN".split()},
}
XL = (0.001, 30)

plt.rcParams.update(
    {
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "pdf.fonttype": 42,
    }
)
fig, axs = plt.subplots(
    1, 4, figsize=(7.4, 5.3), sharey=True, gridspec_kw=dict(wspace=0.08)
)
y = np.arange(len(ORDER))[::-1]
Y = dict(zip(ORDER, y))
for ax, (o, title) in zip(axs, OUTC):
    d = t8[t8.Outcome == o]
    ax.axvspan(XL[0], 1, color="#E64B35", alpha=0.05, lw=0)
    ax.axvline(1, color="black", lw=0.6)
    for sc, lab, mk, col, dy in SC:
        x = d[d.Scenario == sc].set_index("ISO3")
        for iso in ORDER:
            if iso not in x.index or pd.isna(x.loc[iso].iloc[3]):
                continue
            m, lo, hi = x.loc[iso].iloc[3], x.loc[iso].iloc[4], x.loc[iso].iloc[5]
            ax.plot(
                [max(lo, XL[0]), min(hi, XL[1])],
                [Y[iso] + dy] * 2,
                color=col,
                lw=0.7,
                alpha=0.8,
            )
            ax.scatter(
                np.clip(m, XL[0], XL[1]),
                Y[iso] + dy,
                marker=mk,
                s=12,
                color=col,
                edgecolor="white",
                lw=0.3,
                zorder=3,
            )
    if o == "CKD":
        ax.text(
            XL[0] * 1.5,
            Y["NGA"],
            "Not estimable",
            fontsize=5.4,
            color="grey",
            va="center",
        )
    ax.set_xscale("log")
    ax.set_xlim(*XL)
    ax.set_xticks([0.01, 0.1, 1, 10])
    ax.set_xticklabels(["0.01", "0.1", "1", "10"])
    ax.minorticks_off()
    ax.set_title(title, fontsize=7.2, fontweight="bold", pad=4)
    ax.axhline(len(ORDER) - 1.5, color="grey", lw=0.5, ls=(0, (2, 2)))
    ax.grid(axis="y", color="#f2f2f2", lw=0.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
axs[0].set_yticks(y)
axs[0].set_yticklabels(["Global" if i == "GLB" else i for i in ORDER], fontsize=6.6)
axs[0].set_ylim(-0.6, len(ORDER) - 0.4)
for t_, iso in zip(axs[0].get_yticklabels(), ORDER):
    t_.set_color(GRP[iso])
    t_.set_fontweight("bold" if iso == "GLB" else "normal")
fig.supxlabel(
    "Avertable share, dietary policy / GLP-1RA treatment at 20% coverage (log scale)",
    fontsize=7,
    y=0.075,
)
fig.text(
    0.5,
    0.955,
    "← GLP-1RA treatment averts a larger share          Dietary policy averts a larger share →",
    ha="center",
    fontsize=6.2,
    color="grey",
)
fig.legend(
    handles=[
        Line2D([], [], marker=mk, ls="-", lw=0.7, color=col, ms=4, label=lab)
        for _, lab, mk, col, _ in SC
    ],
    loc="lower center",
    ncol=3,
    fontsize=6.2,
    frameon=False,
    bbox_to_anchor=(0.5, 0.0),
)
fig.subplots_adjust(left=0.07, right=0.985, top=0.9, bottom=0.15)
os.makedirs(os.path.dirname(R(a.out)) or ".", exist_ok=True)
fig.savefig(R(a.out) + ".pdf")
fig.savefig(R(a.out) + ".png", dpi=400)
print("saved:", R(a.out))
