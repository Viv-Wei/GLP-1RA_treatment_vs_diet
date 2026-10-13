import subprocess, os, itertools, tempfile, pandas as pd, io

GRID = {
    "policy_maintain": [1.0, 0.8, 0.6, 0.4],
    "wholegrain_disc": [1.0],
    "glp1_retain": [0.37, 0.57, 0.76, 1.13],
    "ramp_diet": [8, 12, 16],
}

from multiprocessing import Pool


def run_one(args):
    pm, wg, rt, rd = args
    rows = []
    env = os.environ | dict(
        SENS_POLICY_MAINTAIN=str(pm), SENS_WHOLEGRAIN=str(wg), SENS_RETAIN=str(rt)
    )
    out = os.path.join(tempfile.gettempdir(), f"sens_{pm}_{wg}_{rt}_{rd}.csv")
    subprocess.run(
        [
            "python",
            "code/05_projection/lifetable_full.py",
            "--years",
            "30",
            "--ramp-diet",
            str(rd),
            "--out",
            out,
        ],
        env=env,
        capture_output=True,
    )
    d = pd.read_csv(out)
    yr = [c for c in d.columns if c.startswith("yr") and c[2:].isdigit()]
    for iso in sorted(d.iso3.unique()):
        sub = d[d.iso3 == iso]
        g = sub.groupby("arm")[yr].sum().T
        g.index = [int(c[2:]) for c in g.index]
        g = g.sort_index()
        if "diet" not in g.columns or "glp1" not in g.columns:
            continue
        cr = g[g.diet > g.glp1]
        tot = g.sum() / 1000
        rows.append(
            dict(
                iso3=iso,
                policy_maintain=pm,
                glp1_retain=rt,
                ramp_diet=rd,
                diet_30yr_k=round(tot.get("diet", 0), 1),
                glp1_30yr_k=round(tot.get("glp1", 0), 1),
                crossover=int(cr.index[0]) if len(cr) else -1,
                diet_wins_30yr=bool(tot.get("diet", 0) > tot.get("glp1", 0)),
            )
        )
    print(f"pm={pm} retain={rt} ramp={rd} done", flush=True)
    return rows


if __name__ == "__main__":
    with Pool(48) as p:
        rows = [
            r
            for part in p.map(run_one, list(itertools.product(*GRID.values())))
            for r in part
        ]
    t = pd.DataFrame(rows)
    t.to_csv("results/sensitivity_sweep_19c.csv", index=False)
    n_sc = t.groupby("iso3").size().iloc[0]
    print(f"\n{t.iso3.nunique()} units x {n_sc} scenarios = {len(t)} rows")
    byc = t.groupby("iso3").diet_wins_30yr.mean().sort_values(ascending=False)
    print("\n=== Share of diet_wins by unit ===")
    print((byc * 100).round(0).astype(int).to_string())
    print(
        f"\nAll-diet units: {(byc==1).sum()}  All-GLP-1 units: {(byc==0).sum()}  Mixed: {((byc>0)&(byc<1)).sum()}"
    )
    cr = t[t.crossover > 0]
    print(
        f"crossover range: {cr.crossover.min()}-{cr.crossover.max()} years, median {cr.crossover.median():.0f} years"
    )
