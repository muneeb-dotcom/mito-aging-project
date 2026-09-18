import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

de = pd.read_csv(r"C:\Users\CC\regen-convergence-v2\results\de_tables\skin_DE_ranked_all_genes.csv")
de = de.rename(columns={"Unnamed: 0": "gene"})
de["gene"] = de["gene"].str.upper()
de = de.groupby("gene")["log2FoldChange"].mean()  # collapse any duplicate symbols

def zscore_long(long_df):
    long_df = long_df.copy()
    long_df["log_tpm"] = np.log2(long_df["tpm"] + 1)
    long_df["z"] = long_df.groupby(["tissue", "gene"])["log_tpm"].transform(
        lambda x: (x - x.mean()) / x.std()
    )
    return long_df

regen_long = zscore_long(pd.read_csv("data/raw/gtex/gtex_regen_signature_long.csv"))

genes_in_data = regen_long["gene"].unique()
matched = [g for g in genes_in_data if g in de.index]
unmatched = [g for g in genes_in_data if g not in de.index]
print(f"Genes matched to DE table for sign: {len(matched)}/{len(genes_in_data)}")
print(f"Unmatched (dropped): {unmatched}")

sign_map = np.sign(de.loc[matched]).to_dict()
n_up = sum(v > 0 for v in sign_map.values())
n_down = sum(v < 0 for v in sign_map.values())
n_zero = sum(v == 0 for v in sign_map.values())
print(f"Up-in-regeneration: {n_up}, up-in-scarring (flipped): {n_down}, zero-FC (dropped): {n_zero}")

regen_long = regen_long[regen_long["gene"].isin(matched)].copy()
regen_long = regen_long[regen_long["gene"].map(sign_map) != 0]
regen_long["signed_z"] = regen_long["z"] * regen_long["gene"].map(sign_map)

regen_scores = regen_long.groupby(["sample_id", "tissue"])["signed_z"].mean().rename("regen_score_v3").reset_index()

mito_scores = pd.read_csv("results/phase2_normalized.csv")[["sample_id", "tissue", "age_mid", "rin", "ischemic_time", "mitonuclear_score_v2"]]
merged = mito_scores.merge(regen_scores, on=["sample_id", "tissue"])
merged.to_csv("results/phase2_signcorrected.csv", index=False)

print("\n--- Sign-corrected regen score vs mitonuclear score: raw Spearman ---")
for tissue, sub in merged.groupby("tissue"):
    rho, p = spearmanr(sub["mitonuclear_score_v2"], sub["regen_score_v3"])
    print(f"{tissue}: n={len(sub)}, rho={rho:.3f}, p={p:.4g}")

print("\n--- Partial (controlling RIN + ischemic time) ---")
for tissue, sub in merged.groupby("tissue"):
    sub = sub.dropna(subset=["rin", "ischemic_time", "mitonuclear_score_v2", "regen_score_v3"])
    X = sm.add_constant(sub[["rin", "ischemic_time"]])
    r1 = sm.OLS(sub["mitonuclear_score_v2"], X).fit().resid
    r2 = sm.OLS(sub["regen_score_v3"], X).fit().resid
    rho, p = spearmanr(r1, r2)
    print(f"{tissue}: partial rho={rho:.3f}, p={p:.4g}")