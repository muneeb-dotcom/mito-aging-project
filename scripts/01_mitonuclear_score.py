import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns

nuclear_oxphos = pd.read_csv("data/raw/mitocarta/nuclear_oxphos_genes.csv")["gene_symbol"].tolist()
mtdna_genes = ["MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
               "MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ATP6","MT-ATP8"]

long = pd.read_csv("data/raw/gtex/gtex_oxphos_long.csv")
print("Loaded:", long.shape)
print("Unique genes in data:", long['gene'].nunique())

wide = long.pivot_table(index=["sample_id", "tissue", "age_mid"], columns="gene", values="tpm").reset_index()

nuclear_present = [g for g in nuclear_oxphos if g in wide.columns]
mtdna_present = [g for g in mtdna_genes if g in wide.columns]
print(f"Nuclear genes matched: {len(nuclear_present)}/{len(nuclear_oxphos)}")
print(f"mtDNA genes matched: {len(mtdna_present)}/{len(mtdna_genes)}")

wide["mitonuclear_score"] = wide[nuclear_present].mean(axis=1) / wide[mtdna_present].mean(axis=1)

wide.to_csv("results/mitonuclear_scores.csv", index=False)

print("\n--- Age correlation per tissue ---")
for tissue, sub in wide.groupby("tissue"):
    rho, p = spearmanr(sub["age_mid"], sub["mitonuclear_score"])
    print(f"{tissue}: n={len(sub)}, Spearman rho={rho:.3f}, p={p:.4g}")

g = sns.lmplot(data=wide, x="age_mid", y="mitonuclear_score", col="tissue",
                col_wrap=3, scatter_kws={"alpha": 0.3})
g.set_axis_labels("Age (bracket midpoint)", "Mitonuclear score (nuclear:mtDNA OXPHOS)")
plt.savefig("figures/mitonuclear_score_vs_age.png", dpi=200, bbox_inches="tight")
print("\nSaved plot to figures/mitonuclear_score_vs_age.png")