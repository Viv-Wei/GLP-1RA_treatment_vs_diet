import os


import sys, pandas as pd

src = sys.argv[1] if len(sys.argv) > 1 else "results/paf_check_vs_official_v2.csv"
dst = sys.argv[2] if len(sys.argv) > 2 else "params/paf_calibration_k.csv"
t = pd.read_csv(src)
t["k"] = (t.official / t.ours).clip(upper=1.5)
k = t.pivot(index="iso3", columns="outcome", values="k").round(4)
k.to_csv(dst)
print(f"=== k = official/ours  ({src} -> {dst}) ===")
print(k.to_string())
print("\nOutcome means:", k.mean().round(3).to_dict())
