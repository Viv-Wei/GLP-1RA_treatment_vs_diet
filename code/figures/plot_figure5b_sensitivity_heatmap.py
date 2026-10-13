import argparse, os
import numpy as np, pandas as pd, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

p = argparse.ArgumentParser()
p.add_argument("--root", default=".")
p.add_argument("--out", default="figures/suppfig_sensitivity_heatmap")
a = p.parse_args()
R = lambda f: os.path.join(a.root, f)

t10 = pd.read_csv(R("results/itt_1012/tableS10_ranking_threshold_itt.csv")).set_index(
    "iso3"
)
sd = pd.read_csv(R("results/itt_1012/sens_deterministic_itt.csv"))
b30 = pd.read_csv(R("results/itt_1012/sens_bmi30_itt.csv")).set_index("iso3")
dc = pd.read_csv(R("results/itt_1012/sens_discount_itt.csv"))
sw = pd.read_csv(R("results/sensitivity_sweep_19c.csv"))
sw["ratio"] = sw.diet_30yr_k / sw.glp1_30yr_k
po = pd.read_csv(
    R("tables/tableS_population_sensitivity_v2.csv"), skiprows=3
).set_index("Country")
fs = R("results/itt_1012/sens_structural_countries_itt.csv")
if os.path.exists(fs):
    stc = pd.read_csv(fs)
    struct = lambda scn: stc[stc.scenario == scn].set_index("iso3").ratio
else:
    stg = pd.read_csv(
        R("tables/tableS_lifetable_structural_sensitivity.csv"), skiprows=3
    ).set_index("Scenario")
    struct = lambda scn: pd.Series({"GLB": stg.loc[scn].iloc[3]})

det = (
    lambda scn, outc="all": sd[(sd.scenario == scn) & (sd.outcomes == outc)]
    .set_index("iso3")
    .ratio
)
sweep = (
    lambda pm, rt, rd: sw[
        np.isclose(sw.policy_maintain, pm)
        & np.isclose(sw.glp1_retain, rt)
        & (sw.ramp_diet == rd)
    ]
    .set_index("iso3")
    .ratio
)
disc = lambda rate: dc[np.isclose(dc.rate, rate)].set_index("iso3").ratio
base = sd[(sd.scenario == "base_0.20") & (sd.outcomes == "all")].set_index("iso3")

ROWS = [
    ("Main analysis", "Main analysis (Monte Carlo median)", t10.threshold_med / 20),
    ("Model structure", "Incidence share of treatment effect 0%", det("split_0")),
    ("Model structure", "Incidence share of treatment effect 100%", det("split_1.0")),
    (
        "Model structure",
        "Same-cycle case fatality for 50% of new cases",
        struct("Same-cycle CFR for 50% of new cases"),
    ),
    (
        "Model structure",
        "Background mortality in diseased ×1.5",
        struct("Background mortality 1.5x for diseased"),
    ),
    (
        "Model structure",
        "Background mortality in diseased ×2.0",
        struct("Background mortality 2.0x for diseased"),
    ),
    (
        "Model structure",
        "Background mortality in diseased ×3.0",
        struct("Background mortality 3.0x for diseased"),
    ),
    (
        "Treatment effect and eligibility",
        "Cardiovascular HR 0.74 (SUSTAIN-6)",
        det("cvhr_0.74"),
    ),
    (
        "Treatment effect and eligibility",
        "Cardiovascular HR 0.86 (SOUL)",
        det("cvhr_0.86"),
    ),
    (
        "Treatment effect and eligibility",
        "CKD HR 0.66 (FLOW lower 95% CI)",
        det("ckdhr_0.66"),
    ),
    (
        "Treatment effect and eligibility",
        "CKD HR 0.88 (FLOW upper 95% CI)",
        det("ckdhr_0.88"),
    ),
    (
        "Treatment effect and eligibility",
        "CKD excluded (IHD and IS only)",
        det("base_0.20", "noCKD"),
    ),
    (
        "Treatment effect and eligibility",
        "Eligibility restricted to BMI ≥30",
        b30.ratio,
    ),
    ("Treatment persistence", "Retention 0.37 (real-world 0.33)", sweep(1.0, 0.37, 8)),
    ("Treatment persistence", "Retention 0.57 (real-world 0.50)", sweep(1.0, 0.57, 8)),
    ("Treatment persistence", "Retention 1.13 (full persistence)", sweep(1.0, 1.13, 8)),
    (
        "Treatment persistence",
        "Retention 0.667 (no ITT correction)",
        base.diet / (base.glp1 * 0.667 / 0.76),
    ),
    (
        "Dietary policy implementation",
        "Policy effect decays to 80%",
        sweep(0.8, 0.76, 8),
    ),
    (
        "Dietary policy implementation",
        "Policy effect decays to 60%",
        sweep(0.6, 0.76, 8),
    ),
    (
        "Dietary policy implementation",
        "Policy effect decays to 40%",
        sweep(0.4, 0.76, 8),
    ),
    ("Dietary policy implementation", "Ramp-up over 12 years", sweep(1.0, 0.76, 12)),
    ("Dietary policy implementation", "Ramp-up over 16 years", sweep(1.0, 0.76, 16)),
    (
        "Dietary policy implementation",
        "Decay 40%, 16-year ramp-up, retention 1.13",
        sweep(0.4, 1.13, 16),
    ),
    ("Population and time", "2053 population structure", po["Ratio 2053"]),
    ("Population and time", "Discounting 3% a year", disc(0.03)),
    ("Population and time", "Discounting 5% a year", disc(0.05)),
]
COLS = ["GLB"] + t10.drop("GLB").sort_values("threshold_med").index.tolist()
M = pd.DataFrame({lab: s.reindex(COLS) for _, lab, s in ROWS}).T
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
M.round(3).to_csv(a.out + "_data.csv")

GCOL = {
    "Main analysis": "#1a1a1a",
    "Model structure": "#7F7F7F",
    "Treatment effect and eligibility": "#E64B35",
    "Treatment persistence": "#B5651D",
    "Dietary policy implementation": "#00A087",
    "Population and time": "#3C5488",
}
GRP_COUNTRY = {
    **{k: "#00A087" for k in "BGR ROU IND MEX ZAF CHN TUR BRA NZL CAN AUS ISR".split()},
    **{k: "#3C5488" for k in "EGY BGD NGA".split()},
    **{k: "#E64B35" for k in "USA DEU GBR JPN".split()},
    "GLB": "#1a1a1a",
}
plt.rcParams.update(
    {
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.5,
        "pdf.fonttype": 42,
    }
)
cmap = LinearSegmentedColormap.from_list(
    "rg", ["#B2182B", "#E64B35", "#F4A582", "#F7F7F7", "#8FD3C1", "#00A087", "#00573F"]
)
norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

HDR = 0.7
FS_LAB, FS_CELL, FS_HDR = 7.0, 5.6, 7.4
rows, y = [], 0
for grp in dict.fromkeys(g for g, _, _ in ROWS):
    if grp != "Main analysis":
        rows.append(("header", grp, y))
        y += HDR
    for g, lab, _ in ROWS:
        if g == grp:
            rows.append(("row", lab, y))
            y += 1
H, W = y, len(COLS)
n_rows = sum(1 for k, _, _ in rows if k == "row")
n_hdr = len(rows) - n_rows
fig = plt.figure(figsize=(7.4, 0.2 * (n_rows + n_hdr) + 1.3))
ax = fig.add_axes([0.335, 0.080, 0.655, 0.856])
for kind, lab, yy in rows:
    if kind == "header":
        ax.text(
            -0.3,
            yy + HDR / 2 + 0.05,
            lab,
            ha="right",
            va="center",
            fontsize=FS_HDR,
            fontweight="bold",
            color=GCOL[lab],
            transform=ax.transData,
            clip_on=False,
        )
        continue
    grp = [g for g, l, _ in ROWS if l == lab][0]
    ax.text(
        -0.3,
        yy + 0.5,
        lab,
        ha="right",
        va="center",
        fontsize=FS_LAB,
        clip_on=False,
        fontweight="bold" if grp == "Main analysis" else "normal",
    )
    for j, iso in enumerate(COLS):
        v = M.loc[lab, iso]
        if pd.isna(v):
            ax.add_patch(
                plt.Rectangle(
                    (j, yy),
                    1,
                    1,
                    facecolor="white",
                    edgecolor="#dddddd",
                    lw=0.3,
                    hatch="////",
                )
            )
            continue
        lv = np.clip(np.log10(v), -1, 1)
        ax.add_patch(
            plt.Rectangle(
                (j, yy), 1, 1, facecolor=cmap(norm(lv)), edgecolor="white", lw=0.4
            )
        )
        txt = f"{v:.2f}" if v < 10 else f"{v:.0f}"
        ax.text(
            j + 0.5,
            yy + 0.52,
            txt,
            ha="center",
            va="center",
            fontsize=FS_CELL,
            color="white" if abs(lv) > 0.55 else "black",
        )
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.set_xticks(np.arange(W) + 0.5)
ax.xaxis.tick_top()
ax.set_xticklabels(
    ["Global" if c == "GLB" else c for c in COLS], rotation=90, fontsize=FS_LAB
)
ax.tick_params(axis="x", pad=1.5)
for t, c in zip(ax.get_xticklabels(), COLS):
    t.set_color(GRP_COUNTRY[c])
    t.set_fontweight("bold" if c == "GLB" else "normal")
ax.axvline(1, color="black", lw=0.8)
ax.set_yticks([])
[s.set_visible(False) for s in ax.spines.values()]
ax.tick_params(length=0)

cax = fig.add_axes([0.47, 0.040, 0.36, 0.011])
cb = matplotlib.colorbar.ColorbarBase(
    cax, cmap=cmap, norm=norm, orientation="horizontal"
)
tk = [0.1, 0.3, 0.5, 1, 2, 3, 10]
cb.set_ticks(np.log10(tk))
cb.set_ticklabels([f"{t:g}" for t in tk])
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=6.5, length=2, pad=1)
cb.set_label(
    "Ratio of deaths averted, dietary policy / GLP-1RA treatment, at 20% coverage",
    fontsize=6.8,
    labelpad=1.5,
)
fig.text(
    0.46,
    0.0455,
    "GLP-1RA treatment\naverts more",
    ha="right",
    va="center",
    fontsize=6.5,
    color="#B2182B",
)
fig.text(
    0.84,
    0.0455,
    "Dietary policy\naverts more",
    ha="left",
    va="center",
    fontsize=6.5,
    color="#00573F",
)
fig.savefig(a.out + ".pdf")
fig.savefig(a.out + ".png", dpi=400)
print("saved:", a.out, "| missing cells:", int(M.isna().sum().sum()))
