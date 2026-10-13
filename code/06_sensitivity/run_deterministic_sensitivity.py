import os, subprocess, tempfile, pandas as pd
from multiprocessing import Pool

SENS_DIR = os.path.join(tempfile.gettempdir(), "sens_itt")
os.makedirs(SENS_DIR, exist_ok=True)
JOBS = [("base_0.20", {}, "0.20"), ("base_0.50", {}, "0.50")]
JOBS += [
    (f"split_{s}", {"SENS_SPLIT_INC": str(s)}, "0.20")
    for s in [0, 0.1, 0.2, 0.33, 0.5, 0.67, 1.0]
]
JOBS += [
    (f"cvhr_{h}", {"SENS_CV_HR": str(h)}, "0.20")
    for h in [0.74, 0.76, 0.78, 0.82, 0.84, 0.86]
]
JOBS += [(f"ckdhr_{h}", {"SENS_CKD_HR": str(h)}, "0.20") for h in [0.66, 0.88]]


def run(job):
    name, extra, cov = job
    out = os.path.join(SENS_DIR, f"{name}.csv")
    r = subprocess.run(
        [
            "python",
            "code/05_projection/lifetable_full.py",
            "--years",
            "30",
            "--coverage",
            cov,
            "--diet-scenario",
            "feasible",
            "--out",
            out,
        ],
        env=os.environ | extra,
        capture_output=True,
        text=True,
    )
    if r.returncode:
        return name, None
    t = pd.read_csv(out)
    yr = [c for c in t.columns if c.startswith("yr")]
    t["cum30"] = t[yr].sum(axis=1)
    res = []
    for lab, keep in [
        ("all", ["IHD", "IschStroke", "CKD"]),
        ("noCKD", ["IHD", "IschStroke"]),
    ]:
        s = (
            t[t.outcome.isin(keep)]
            .groupby(["iso3", "arm"])
            .cum30.sum()
            .unstack("arm")
            .reset_index()
        )
        s["scenario"], s["outcomes"], s["coverage"] = name, lab, cov
        res.append(s)
    return name, pd.concat(res)


if __name__ == "__main__":
    with Pool(len(JOBS)) as p:
        out = p.map(run, JOBS)
    bad = [n for n, d in out if d is None]
    if bad:
        print("failed: ", bad)
    r = pd.concat(d for n, d in out if d is not None)
    r["ratio"] = r.diet / r.glp1
    r["add_diet_pct"] = (r.combo - r.glp1) / r.glp1 * 100
    r.to_csv("results/itt_1012/sens_deterministic_itt.csv", index=False)
    a = r[r.outcomes == "all"]
    summ = a.groupby("scenario").apply(
        lambda x: pd.Series(
            dict(
                glb_ratio=x[x.iso3 == "GLB"].ratio.iloc[0],
                glb_threshold=20 * x[x.iso3 == "GLB"].ratio.iloc[0],
                diet_leads=(x[x.iso3 != "GLB"].ratio > 1).sum(),
            )
        ),
        include_groups=False,
    )
    base = a[a.scenario == "base_0.20"].set_index("iso3").ratio > 1
    summ["flipped_vs_base"] = [
        ", ".join(
            sorted(
                (a[a.scenario == s].set_index("iso3").ratio > 1)
                .pipe(lambda x: x[x != base])
                .index
            )
        )
        or "—"
        for s in summ.index
    ]
    pd.set_option("display.width", 200)
    print(summ.round(2).to_string())
    nc = r[(r.outcomes == "noCKD") & r.scenario.str.startswith("base")]
    print(
        "\nExcluding CKD:",
        nc.groupby("scenario")
        .apply(
            lambda x: f"global ratio {x[x.iso3=='GLB'].ratio.iloc[0]:.2f}, diet leads {(x[x.iso3!='GLB'].ratio>1).sum()}/19",
            include_groups=False,
        )
        .to_dict(),
    )
