import argparse, os
import numpy as np, pandas as pd, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

p = argparse.ArgumentParser()
p.add_argument("--root", default=".")
p.add_argument("--xlsx", default="tables/Supplementary_Tables_final1012_OW.xlsx")
p.add_argument("--out", default="figures/suppfig2_mediation")
a = p.parse_args()
R = lambda f: os.path.join(a.root, f)
tab = lambda n: pd.read_excel(R(a.xlsx), f"Supplementary Table {n}", header=1).iloc[:-1]

t3, t4, t10 = tab(3), tab(4), tab(10)
name2iso = dict(zip(t4.Country, t4.ISO3))
t10["iso3"] = t10.Country.map(name2iso)
ORDER = ["GLB"] + t10[t10.iso3 != "GLB"].sort_values(t10.columns[6]).iso3.tolist()
OUTC = [
    ("Ischaemic heart disease", "IHD"),
    ("Ischaemic stroke", "IS"),
    ("Intracerebral haemorrhage", "ICH"),
    ("Subarachnoid haemorrhage", "SAH"),
    ("Type 2 diabetes", "T2D"),
    ("Chronic kidney disease", "CKD"),
]
CATS = [
    ("Fully mediated (%)", "Fully mediated (>99%)", "#00A087"),
    ("Partially mediated (%)", "Partially mediated (10–50%)", "#91D1C2"),
    ("Minimally mediated (%)", "Minimally mediated (<10%)", "#F39B7F"),
    ("Unmediated (%)", "Unmediated (no quantified pathway)", "#E64B35"),
]
GRP = {
    "GLB": "#1a1a1a",
    **{k: "#00A087" for k in "BGR ROU IND MEX ZAF CHN TUR BRA NZL CAN AUS ISR".split()},
    **{k: "#3C5488" for k in "EGY BGD NGA".split()},
    **{k: "#E64B35" for k in "USA DEU GBR JPN".split()},
}
assert set(o for o, _ in OUTC) == set(t3.Outcome.unique()), t3.Outcome.unique()

plt.rcParams.update(
    {
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "pdf.fonttype": 42,
    }
)
fig = plt.figure(figsize=(7.4, 5.4))
y = np.arange(len(ORDER))[::-1]
Y = dict(zip(ORDER, y))

L, W, G = 0.065, 0.088, 0.008
axs = [fig.add_axes([L + i * (W + G), 0.17, W, 0.76]) for i in range(len(OUTC))]
for ax, (outc, short) in zip(axs, OUTC):
    d = t3[t3.Outcome == outc].set_index("ISO3")
    left = np.zeros(len(ORDER))
    for col, lab, c in CATS:
        v = np.array(
            [float(d.loc[iso, col]) if iso in d.index else 0.0 for iso in ORDER]
        )
        ax.barh(y, v, left=left, height=0.78, color=c, edgecolor="white", lw=0.25)
        left += v
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100] if ax is axs[-1] else [0, 50])
    ax.tick_params(axis="x", labelsize=5.8)
    ax.set_ylim(-0.6, len(ORDER) - 0.4)
    ax.set_yticks(y)
    ax.tick_params(axis="y", length=0)
    ax.set_title(short, fontsize=7.2, fontweight="bold", pad=3)
    ax.axhline(len(ORDER) - 1.5, color="white", lw=1.6)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_yticklabels(
        ["Global" if i == "GLB" else i for i in ORDER] if ax is axs[0] else [],
        fontsize=6.4,
    )
for t_, iso in zip(axs[0].get_yticklabels(), ORDER):
    t_.set_color(GRP[iso])
    t_.set_fontweight("bold" if iso == "GLB" else "normal")
fig.text(
    L + 3 * (W + G) - G / 2,
    0.095,
    "Share of diet-attributable burden by mediation category (%)",
    ha="center",
    fontsize=6.8,
)
fig.legend(
    handles=[Patch(color=c, label=l) for _, l, c in CATS],
    loc="lower center",
    bbox_to_anchor=(0.34, -0.005),
    ncol=2,
    fontsize=6.0,
    frameon=False,
    handlelength=1.1,
    columnspacing=1.0,
)

axB = fig.add_axes([0.70, 0.17, 0.255, 0.76])
RULES = [
    ("Lenient rule (%)", "Lenient", "o", "#0F6E56", 34),
    (
        "Weighted rule, main analysis (%)",
        "Weighted (main analysis)",
        "s",
        "#1D9E75",
        30,
    ),
    ("Strict rule (%)", "Strict", "D", "#9FE1CB", 26),
]
d4 = t4.set_index("ISO3")
for iso in ORDER:
    vals = [float(d4.loc[iso, c]) for c, *_ in RULES]
    axB.plot([min(vals), max(vals)], [Y[iso]] * 2, color="#888780", lw=0.8, zorder=1)
    for (c, lab, mk, col, s), v in zip(RULES, vals):
        axB.scatter(
            v, Y[iso], marker=mk, s=s, color=col, edgecolor="#0F6E56", lw=0.5, zorder=3
        )
axB.set_xlim(20, 100)
axB.set_xticks([20, 40, 60, 80, 100])
axB.set_ylim(-0.6, len(ORDER) - 0.4)
axB.set_yticks(y)
axB.set_yticklabels([])
axB.tick_params(axis="y", length=0)
axB.axhline(len(ORDER) - 1.5, color="grey", lw=0.5, ls=(0, (2, 2)))
axB.grid(axis="x", color="#eeeeee", lw=0.5)
axB.set_axisbelow(True)
for s in ["top", "right"]:
    axB.spines[s].set_visible(False)
axB.set_xlabel(
    "Share of diet-attributable cardiovascular deaths\noutside SBP, FPG, and LDL-C (%)",
    fontsize=6.8,
)
axB.legend(
    handles=[
        Line2D(
            [], [], marker=mk, ls="", mfc=col, mec="#0F6E56", ms=np.sqrt(s), label=lab
        )
        for _, lab, mk, col, s in RULES
    ],
    loc="upper center",
    bbox_to_anchor=(0.42, -0.13),
    ncol=3,
    columnspacing=0.6,
    handletextpad=0.1,
    fontsize=5.8,
    frameon=False,
    title="Accounting rule",
    title_fontsize=6.0,
)
fig.text(0.015, 0.95, "a", fontsize=12, fontweight="bold")
fig.text(0.665, 0.95, "b", fontsize=12, fontweight="bold")
os.makedirs(os.path.dirname(R(a.out)) or ".", exist_ok=True)
fig.savefig(R(a.out) + ".pdf")
fig.savefig(R(a.out) + ".png", dpi=400)
print("saved:", R(a.out), "| global weighted rule:", float(d4.loc["GLB", RULES[1][0]]))
