import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns

mito = pd.read_csv("results/mitonuclear_scores.csv")[["sample_id", "tissue", "age_mid", "mitonuclear_score"]]

regen_long = pd.read_csv("data/raw/gtex/gtex_regen_signature_long.csv")
matched_genes = sorted(regen_long["gene"].unique())
print(f"Regen signature genes matched in GTEx: {len(matched_genes)}")

regen_wide = regen_long.pivot_table(index="sample_id", columns="gene", values="tpm")
regen_wide["regen_score"] = regen_wide.mean(axis=1)
regen_wide = regen_wide[["regen_score"]].reset_index()

merged = mito.merge(regen_wide, on="sample_id", how="inner")
merged.to_csv("results/phase2_merged.csv", index=False)
print(f"Merged samples: {len(merged)}")

print("\n--- Mitonuclear score vs regen signature, per tissue ---")
for tissue, sub in merged.groupby("tissue"):
    rho, p = spearmanr(sub["mitonuclear_score"], sub["regen_score"])
    print(f"{tissue}: n={len(sub)}, Spearman rho={rho:.3f}, p={p:.4g}")

print("\n--- Mixed-effects model (pooled, tissue as random effect) ---")
model = smf.mixedlm("regen_score ~ mitonuclear_score", data=merged, groups=merged["tissue"])
result = model.fit()
print(result.summary())

g = sns.lmplot(data=merged, x="mitonuclear_score", y="regen_score", col="tissue",
                col_wrap=3, scatter_kws={"alpha": 0.3})
g.set_axis_labels("Mitonuclear score", "Regeneration signature score")
plt.savefig("figures/phase2_mitonuclear_vs_regen.png", dpi=200, bbox_inches="tight")
print("\nSaved plot to figures/phase2_mitonuclear_vs_regen.png")