import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

def zscore_long(long_df):
    long_df = long_df.copy()
    long_df["log_tpm"] = np.log2(long_df["tpm"] + 1)
    long_df["z"] = long_df.groupby(["tissue", "gene"])["log_tpm"].transform(
        lambda x: (x - x.mean()) / x.std()
    )
    return long_df

mtdna_genes = ["MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
               "MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ATP6","MT-ATP8"]

mito_long = zscore_long(pd.read_csv("data/raw/gtex/gtex_oxphos_long.csv"))
nuclear_z = mito_long[~mito_long["gene"].isin(mtdna_genes)].groupby(["sample_id","tissue"])["z"].mean().rename("nuclear_z")
mtdna_z = mito_long[mito_long["gene"].isin(mtdna_genes)].groupby(["sample_id","tissue"])["z"].mean().rename("mtdna_z")

mito_scores = pd.concat([nuclear_z, mtdna_z], axis=1).reset_index()
mito_scores["mitonuclear_score_v2"] = mito_scores["nuclear_z"] - mito_scores["mtdna_z"]

regen_long = zscore_long(pd.read_csv("data/raw/gtex/gtex_regen_signature_long.csv"))
regen_scores = regen_long.groupby(["sample_id","tissue"])["z"].mean().rename("regen_score_v2").reset_index()

old = pd.read_csv("results/phase2_with_qc.csv")[["sample_id","tissue","age_mid","rin","ischemic_time"]]
merged = old.merge(mito_scores, on=["sample_id","tissue"]).merge(regen_scores, on=["sample_id","tissue"])
merged.to_csv("results/phase2_normalized.csv", index=False)

print("--- Normalized scores: raw Spearman, per tissue ---")
for tissue, sub in merged.groupby("tissue"):
    rho, p = spearmanr(sub["mitonuclear_score_v2"], sub["regen_score_v2"])
    print(f"{tissue}: n={len(sub)}, rho={rho:.3f}, p={p:.4g}")

print("\n--- Normalized scores: partial (controlling RIN + ischemic time) ---")
for tissue, sub in merged.groupby("tissue"):
    sub = sub.dropna(subset=["rin","ischemic_time","mitonuclear_score_v2","regen_score_v2"])
    X = sm.add_constant(sub[["rin","ischemic_time"]])
    r1 = sm.OLS(sub["mitonuclear_score_v2"], X).fit().resid
    r2 = sm.OLS(sub["regen_score_v2"], X).fit().resid
    rho, p = spearmanr(r1, r2)
    print(f"{tissue}: partial rho={rho:.3f}, p={p:.4g}")