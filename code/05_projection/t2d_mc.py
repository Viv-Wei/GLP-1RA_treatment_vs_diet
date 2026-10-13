import argparse, os, math
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--draws", type=int, default=500)
ap.add_argument("--years", type=int, default=30)
ap.add_argument("--ramp-glp1", type=int, default=1)
ap.add_argument("--ramp-diet", type=int, default=8)
ap.add_argument("--diet-scenario", default="feasible")
A = ap.parse_args()

cm = pd.read_csv("params/country_map.csv")
ISO = dict(zip(cm.gbd_location_name, cm.iso3))
AGE15 = [f"{x}-{x+4} years" for x in range(25, 95, 5)] + ["95+ years"]
ramp = lambda t, yrs: 1.0 / (1.0 + math.exp(-(4.0 / yrs) * (t - yrs / 1.6)))
d = pd.read_csv("data/gbd/cause_totals/queryB_t2d.csv", low_memory=False)
d = d[
    (d.measure_name == "Incidence")
    & (d.metric_name == "Number")
    & (d.year == 2023)
    & d.age_name.isin(AGE15)
    & d.sex_name.isin(["Male", "Female"])
    & d.location_name.isin(ISO)
]
inc0 = d.assign(iso3=d.location_name.map(ISO)).groupby("iso3").val.sum()
ow = pd.read_csv("params/overweight_prevalence.csv").set_index("iso3").overweight_prev
dp = pd.read_csv(os.environ.get("LT_DP_OVERRIDE", "results/diet_policy_arm.csv"))
diet_pif = (
    dp[(dp.scenario == A.diet_scenario) & (dp.outcome == "T2D")]
    .groupby("iso3")
    .averted_frac_alldeaths.mean()
)
tp = pd.read_csv("params/glp1_trial_params.csv")
r = tp[tp.param == "t2d_incidence_hr"].iloc[0]
MU, SD = np.log(r.value), (np.log(r.ci_hi) - np.log(r.ci_lo)) / (2 * 1.96)
aff = pd.read_csv("params/coverage_afford_cap.csv")
AFF = dict(zip(aff.iloc[:, 0], aff.iloc[:, 1]))

isos = list(inc0.index)
R_G = sum(ramp(t, A.ramp_glp1) for t in range(1, A.years + 1))
R_D = sum(ramp(t, A.ramp_diet) for t in range(1, A.years + 1))
rng = np.random.default_rng(4200)
hr = np.exp(rng.normal(MU, SD, A.draws))
rw = rng.uniform(0.50, 0.83, A.draws)
p_stop = (1 - rw) * 1.5
ret = (1 - 0.75 * p_stop) / 0.869
dfac = rng.normal(1.0, 0.06, (A.draws, len(isos)))

base = inc0.values
owv = np.array([ow.get(i, np.nan) for i in isos])
dpv = np.array([diet_pif.get(i, 0.0) for i in isos])
rows = []
diet = base * dpv * dfac * R_D
rows += [
    dict(
        iso3=i,
        coverage="any",
        arm="Dietary policy (feasible)",
        med=np.median(diet[:, k]),
        lo=np.quantile(diet[:, k], 0.025),
        hi=np.quantile(diet[:, k], 0.975),
    )
    for k, i in enumerate(isos)
]
for cov_lab in ["20% coverage", "50% coverage", "Maximum affordable coverage"]:
    cov = np.array(
        [
            (
                0.20
                if cov_lab.startswith("20")
                else 0.50 if cov_lab.startswith("50") else AFF.get(i, 1.0)
            )
            for i in isos
        ]
    )
    glp = base * cov * owv * ((1 - hr) * ret)[:, None] * R_G
    rows += [
        dict(
            iso3=i,
            coverage=cov_lab,
            arm="GLP-1RA treatment",
            med=np.median(glp[:, k]),
            lo=np.quantile(glp[:, k], 0.025),
            hi=np.quantile(glp[:, k], 0.975),
        )
        for k, i in enumerate(isos)
    ]
    rat = diet / glp
    rows += [
        dict(
            iso3=i,
            coverage=cov_lab,
            arm="Ratio diet/GLP-1RA",
            med=np.median(rat[:, k]),
            lo=np.quantile(rat[:, k], 0.025),
            hi=np.quantile(rat[:, k], 0.975),
        )
        for k, i in enumerate(isos)
    ]
out = pd.DataFrame(rows)
os.makedirs("results/itt_1012", exist_ok=True)
out.to_csv("results/itt_1012/t2d_mc_summary.csv", index=False)
g = out[out.iso3 == "GLB"]
print(
    f"T2D Monte Carlo, {A.draws} draws (HR median {np.median(hr):.3f}, retention median {np.median(ret):.3f})"
)
print(
    (
        g.assign(
            med=g.med.where(g.arm.str.startswith("Ratio"), g.med / 1e6),
            lo=g.lo.where(g.arm.str.startswith("Ratio"), g.lo / 1e6),
            hi=g.hi.where(g.arm.str.startswith("Ratio"), g.hi / 1e6),
        )
    )
    .round(2)
    .to_string(index=False)
)
