import argparse, os
import numpy as np, pandas as pd, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

p = argparse.ArgumentParser()
p.add_argument("--root", default=".")
p.add_argument("--out", default="figures/figure3_threshold")
a = p.parse_args()
R = lambda f: os.path.join(a.root, f)

t9 = pd.read_csv(R("results/itt_1012/tableS9_deaths_itt.csv"))
w = t9.pivot_table(index="iso3", columns=["arm", "coverage"], values="med")
d = pd.DataFrame(
    {
        "diet": w[("Dietary policy (feasible)", "20% coverage")],
        "asp": w[("Dietary policy (aspirational)", "20% coverage")],
        "g20": w[("GLP-1RA treatment", "20% coverage")],
        "gmax": w[("GLP-1RA treatment", "Maximum affordable coverage")],
    }
)
d["maxcov"] = 20 * d.gmax / d.g20
d["thr_asp"] = 20 * d.asp / d.g20

t10 = pd.read_csv(R("results/itt_1012/tableS10_ranking_threshold_itt.csv")).set_index(
    "iso3"
)
d = d.join(
    t10[["threshold_med", "threshold_lo", "threshold_hi"]].set_axis(
        ["thr", "lo", "hi"], axis=1
    )
)

sw = pd.read_csv(R("results/sensitivity_sweep_19c.csv"))
sw = sw[(sw.policy_maintain == 1.0) & (sw.ramp_diet == 8)]
for rt in [0.37, 1.13]:
    x = sw[np.isclose(sw.glp1_retain, rt)].set_index("iso3")
    d[f"thr_r{rt}"] = 20 * x.diet_30yr_k / x.glp1_30yr_k

GRP = {
    "GLB": "Global",
    **{
        k: "Patent expiry ≤2027"
        for k in "BGR ROU IND MEX ZAF CHN TUR BRA NZL CAN AUS ISR".split()
    },
    **{k: "Never filed" for k in "EGY BGD NGA".split()},
    **{k: "Protected" for k in "USA DEU GBR JPN".split()},
}
COL = {
    "Global": "#1a1a1a",
    "Patent expiry ≤2027": "#00A087",
    "Never filed": "#3C5488",
    "Protected": "#E64B35",
}
d["grp"] = d.index.map(GRP)
assert d.notna().all().all(), d[d.isna().any(axis=1)]
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
d.round(2).to_csv(a.out + "_data.csv")

plt.rcParams.update(
    {
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "hatch.linewidth": 0.4,
        "pdf.fonttype": 42,
    }
)
fig = plt.figure(figsize=(7.4, 5.2))
axA = fig.add_axes([0.065, 0.25, 0.46, 0.67])
axB = fig.add_axes([0.62, 0.25, 0.36, 0.67])

c = np.linspace(3, 100, 400)
axA.axhspan(1, 40, color="#00A087", alpha=0.05, lw=0)
axA.axhspan(0.01, 1, color="#E64B35", alpha=0.05, lw=0)
axA.text(
    99,
    30,
    "Dietary policy averts more",
    ha="right",
    va="top",
    color="#00795f",
    fontsize=7,
    style="italic",
)
axA.text(
    99,
    0.0115,
    "GLP-1RA treatment averts more",
    ha="right",
    va="bottom",
    color="#b23a29",
    fontsize=7,
    style="italic",
)
for x0, lab in [(20, "Near-term\n(20%)"), (50, "50%")]:
    axA.axvline(x0, color="grey", lw=0.6, ls=(0, (3, 2)), zorder=1)
    axA.text(x0 + 1, 30, lab, va="top", fontsize=6.5, color="grey")
axA.axhline(1, color="black", lw=0.8, zorder=2)

g = d.loc["GLB"]
axA.fill_between(c, g.lo / c, g.hi / c, color="grey", alpha=0.18, lw=0, zorder=1)
axA.annotate(
    f"Global parity at\n{g.thr:.0f}% ({g.lo:.0f}–{g.hi:.0f}%)",
    xy=(g.thr, 1),
    xytext=(g.thr + 9, 2.6),
    fontsize=6.5,
    fontweight="bold",
    arrowprops=dict(arrowstyle="->", lw=0.6, color="black"),
)

ends = []
for iso, r in d.iterrows():
    col, lw, z = COL[r.grp], (1.8 if iso == "GLB" else 0.9), (5 if iso == "GLB" else 3)
    m = min(r.maxcov, 100)
    cin, cout = c[c <= m + 1e-9], c[c >= m - 1e-9]
    axA.plot(
        cin, r.thr / cin, color=col, lw=lw, alpha=1 if iso == "GLB" else 0.85, zorder=z
    )
    if len(cout) > 1:
        axA.plot(
            cout,
            r.thr / cout,
            color=col,
            lw=lw * 0.8,
            ls=(0, (1, 1.5)),
            alpha=0.6,
            zorder=z,
        )
    if r.thr <= 100:
        axA.scatter(
            r.thr,
            1,
            s=22 if iso == "GLB" else 12,
            color=col,
            edgecolor="white",
            lw=0.5,
            zorder=6,
        )
    ends.append([iso, r.thr / 100, col])

ends = sorted(ends, key=lambda x: x[1])
ys = np.log10([e[1] for e in ends])
for i in range(1, len(ys)):
    ys[i] = max(ys[i], ys[i - 1] + 0.085)
for (iso, y0, col), yl in zip(ends, ys):
    axA.annotate(
        "Global" if iso == "GLB" else iso,
        xy=(100, y0),
        xytext=(103, 10**yl),
        fontsize=6,
        color=col,
        va="center",
        fontweight="bold" if iso == "GLB" else "normal",
        annotation_clip=False,
        arrowprops=dict(
            arrowstyle="-", lw=0.3, color=col, alpha=0.6, shrinkA=0, shrinkB=0
        ),
    )

axA.set_yscale("log")
axA.set_ylim(0.01, 40)
axA.set_xlim(0, 100)
yt = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]
axA.set_yticks(yt)
axA.set_yticklabels([f"{v:g}" for v in yt])
axA.set_xlabel("GLP-1RA treatment coverage (% of adults with BMI ≥25)")
axA.set_ylabel("Deaths averted, dietary policy / GLP-1RA treatment (log scale)")
axA.spines[["top", "right"]].set_visible(False)
axA.text(-0.09, 1.04, "a", transform=axA.transAxes, fontsize=12, fontweight="bold")
axA.legend(
    handles=[Line2D([], [], color=v, lw=1.4, label=k) for k, v in COL.items()],
    title="Semaglutide patent status",
    title_fontsize=6.5,
    fontsize=6.3,
    loc="lower left",
    bbox_to_anchor=(0.0, 0.07),
    frameon=False,
    handlelength=1.6,
)

order = ["GLB"] + d.drop("GLB").sort_values("thr").index.tolist()
y = np.arange(len(order))[::-1]
XM = 150
for x0 in [20, 50]:
    axB.axvline(x0, color="grey", lw=0.6, ls=(0, (3, 2)), zorder=1)
axB.axvspan(0, 20, color="#E64B35", alpha=0.06, lw=0)
for yy, iso in zip(y, order):
    r = d.loc[iso]
    col = COL[r.grp]
    lo_p, hi_p = sorted([r["thr_r1.13"], r["thr_r0.37"]])
    axB.plot(
        [lo_p, min(hi_p, XM)],
        [yy, yy],
        color=col,
        lw=4.2,
        alpha=0.18,
        solid_capstyle="butt",
        zorder=2,
    )
    axB.plot([r.lo, min(r.hi, XM)], [yy, yy], color=col, lw=0.9, zorder=3)
    axB.scatter(
        min(r.thr, XM),
        yy,
        s=26 if iso == "GLB" else 18,
        color=col,
        zorder=5,
        edgecolor="white",
        lw=0.5,
    )
    axB.scatter(
        min(r.thr_asp, XM),
        yy,
        s=14,
        marker="D",
        facecolor="white",
        edgecolor=col,
        lw=0.7,
        zorder=4,
    )
    m = min(r.maxcov, 100)
    axB.add_patch(
        plt.Rectangle(
            (m, yy - 0.42),
            XM + 6 - m,
            0.84,
            facecolor="none",
            edgecolor="#9a9a9a",
            hatch="////",
            lw=0,
            zorder=1.5,
        )
    )
    if max(hi_p, r.hi, r.thr_asp) > XM:
        axB.annotate(
            "",
            xy=(XM + 4, yy),
            xytext=(XM - 1, yy),
            arrowprops=dict(arrowstyle="->", lw=0.6, color=col),
            annotation_clip=False,
        )
axB.set_yticks(y)
axB.set_yticklabels(["Global" if i == "GLB" else i for i in order])
for t, iso in zip(axB.get_yticklabels(), order):
    t.set_color(COL[d.loc[iso, "grp"]])
    t.set_fontweight("bold" if iso == "GLB" else "normal")
axB.set_xlim(0, XM)
axB.set_ylim(-0.7, len(order) - 0.3)
axB.set_xticks([0, 20, 50, 100, 150])
axB.set_xlabel("Coverage needed to match dietary policy\n(% of adults with BMI ≥25)")
axB.spines[["top", "right"]].set_visible(False)
axB.text(-0.2, 1.04, "b", transform=axB.transAxes, fontsize=12, fontweight="bold")
axB.legend(
    handles=[
        Line2D(
            [], [], marker="o", ls="", color="grey", ms=4.5, label="Median (95% UI)"
        ),
        Line2D(
            [],
            [],
            color="grey",
            lw=4.2,
            alpha=0.3,
            label="Range across persistence\n(retention 0.37–1.13)",
        ),
        Line2D(
            [],
            [],
            marker="D",
            ls="",
            mfc="white",
            mec="grey",
            ms=3.6,
            label="Versus aspirational scenario",
        ),
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor="white",
            edgecolor="#9a9a9a",
            hatch="////",
            lw=0,
            label="Beyond maximum affordable\ncoverage at generic prices",
        ),
    ],
    loc="upper center",
    bbox_to_anchor=(0.5, -0.17),
    ncol=2,
    fontsize=6,
    frameon=False,
    handlelength=1.6,
)

fig.savefig(a.out + ".pdf")
fig.savefig(a.out + ".png", dpi=400)
print("saved:", a.out + ".pdf / .png / _data.csv")
