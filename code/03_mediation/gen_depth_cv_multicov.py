import pandas as pd, json, glob, os

os.chdir(os.environ.get("CRA_ROOT", os.getcwd()))
O2S = {"IHD": "IHD", "Ischaemic stroke": "IschStroke", "ICH": "ICH", "SAH": "SAH"}
COV = {
    "cov020": ("params", "20%"),
    "cov050": ("params/depth_cov050", "50%"),
    "afford": ("params/depth_afford", "affordable"),
}
SEG_NP = "Unreachable (no pathway)"
SEG_OPEN = "Pathway open, flow insufficient"

isos = sorted(
    os.path.basename(f)[len("figure3_triple_") : -4]
    for f in glob.glob("results/figure3_triple_*.csv")
    if "depth" not in f
)
country, factor_rows = [], []
for iso in isos:
    t = pd.read_csv(f"results/figure3_triple_{iso}.csv")
    t = t[t.outcome.isin(O2S) & (t.deaths > 0)]
    for tag, (ddir, lab) in COV.items():
        seg_off = f"Offset by GLP-1 ({lab} coverage)"
        rows, missing = [], 0
        for _, r in t.iterrows():
            o = O2S[r.outcome]
            if r.mediator == "Direct":
                rows.append((r.factor, r.mediator, r.outcome, SEG_NP, r.deaths))
                continue
            fp = f"{ddir}/figure3_depth_{iso}_{o}.json"
            if not os.path.exists(fp):
                missing += 1
                dp = 0.0
            else:
                dp = max(float(json.load(open(fp)).get(r.mediator, -1.0)), 0.0)
            dp = min(dp, 1.0)
            if dp > 0:
                rows.append((r.factor, r.mediator, r.outcome, seg_off, r.deaths * dp))
            if dp < 1:
                rows.append(
                    (r.factor, r.mediator, r.outcome, SEG_OPEN, r.deaths * (1 - dp))
                )
        d = pd.DataFrame(
            rows, columns=["factor", "mediator", "outcome", "seg", "deaths"]
        )
        assert abs(d.deaths.sum() - t.deaths.sum()) < max(
            1, 1e-9 * t.deaths.sum()
        ), f"{iso} {tag} not conserved"
        if missing:
            print(
                f"[warn] {iso} {tag}: {missing} flows without depth file -> treated as not offset"
            )
        d.to_csv(f"results/f3cv_flows_{tag}_{iso}.csv", index=False)
        fb = d.groupby(["factor", "seg"]).deaths.sum().reset_index()
        fb.assign(med=fb.deaths, lo=fb.deaths, hi=fb.deaths)[
            ["factor", "seg", "med", "lo", "hi"]
        ].to_csv(f"results/f3cv_depth_{tag}_{iso}.csv", index=False)
        s = d.groupby("seg").deaths.sum()
        tot = s.sum()
        country.append(
            dict(
                iso3=iso,
                coverage=lab,
                total_cv_deaths=round(tot),
                offset_pct=round(s.get(seg_off, 0) / tot * 100, 1),
                mediated_not_offset_pct=round(s.get(SEG_OPEN, 0) / tot * 100, 1),
                no_pathway_pct=round(s.get(SEG_NP, 0) / tot * 100, 1),
            )
        )
        for fac, g in d.groupby("factor"):
            ft = g.deaths.sum()
            factor_rows.append(
                dict(
                    iso3=iso,
                    coverage=lab,
                    factor=fac,
                    attributable_cv_deaths=round(ft),
                    offset_pct=round(g[g.seg == seg_off].deaths.sum() / ft * 100, 1),
                    mediated_not_offset_pct=round(
                        g[g.seg == SEG_OPEN].deaths.sum() / ft * 100, 1
                    ),
                    no_pathway_pct=round(g[g.seg == SEG_NP].deaths.sum() / ft * 100, 1),
                )
            )

C = pd.DataFrame(country)
C.to_csv("results/depth_cv_country.csv", index=False)
F = pd.DataFrame(factor_rows)
F.to_csv("results/depth_cv_factor.csv", index=False)
print("countries:", len(isos))
print("\n=== GLB, CV only ===\n", C[C.iso3 == "GLB"].to_string(index=False))
print("\n=== GLB by factor ===")
print(
    F[F.iso3 == "GLB"]
    .pivot(index="factor", columns="coverage", values="offset_pct")
    .join(
        F[(F.iso3 == "GLB") & (F.coverage == "20%")].set_index("factor")[
            ["attributable_cv_deaths", "no_pathway_pct"]
        ]
    )
    .sort_values("attributable_cv_deaths", ascending=False)
    .to_string()
)
