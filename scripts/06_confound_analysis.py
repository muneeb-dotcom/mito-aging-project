import pandas as pd
from scipy.stats import spearmanr

merged = pd.read_csv("results/phase2_merged.csv")
qc = pd.read_csv("data/raw/gtex/gtex_qc_metadata.csv")

df = merged.merge(qc, on=["sample_id", "tissue"], how="inner")
print(f"Merged with QC: {len(df)} samples")
print(df[["rin", "ischemic_time"]].describe())

print("\n--- Does RIN/ischemic time correlate with our two scores? ---")
for tissue, sub in df.groupby("tissue"):
    r1, p1 = spearmanr(sub["rin"], sub["mitonuclear_score"])
    r2, p2 = spearmanr(sub["rin"], sub["regen_score"])
    r3, p3 = spearmanr(sub["ischemic_time"], sub["mitonuclear_score"])
    r4, p4 = spearmanr(sub["ischemic_time"], sub["regen_score"])
    print(f"\n{tissue} (n={len(sub)}):")
    print(f"  RIN vs mitonuclear_score:      rho={r1:.3f}, p={p1:.4g}")
    print(f"  RIN vs regen_score:            rho={r2:.3f}, p={p2:.4g}")
    print(f"  ischemic_time vs mitonuclear:  rho={r3:.3f}, p={p3:.4g}")
    print(f"  ischemic_time vs regen_score:  rho={r4:.3f}, p={p4:.4g}")

print("\n--- Partial correlation: mitonuclear vs regen, controlling for RIN + ischemic time ---")
import statsmodels.api as sm

for tissue, sub in df.groupby("tissue"):
    sub = sub.dropna(subset=["rin", "ischemic_time", "mitonuclear_score", "regen_score"])
    X_mito = sm.add_constant(sub[["rin", "ischemic_time"]])
    resid_mito = sm.OLS(sub["mitonuclear_score"], X_mito).fit().resid
    resid_regen = sm.OLS(sub["regen_score"], X_mito).fit().resid
    rho, p = spearmanr(resid_mito, resid_regen)
    print(f"{tissue}: partial rho={rho:.3f}, p={p:.4g}")

df.to_csv("results/phase2_with_qc.csv", index=False)