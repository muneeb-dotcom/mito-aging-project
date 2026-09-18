import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import shap

def zscore_long(long_df):
    long_df = long_df.copy()
    long_df["log_tpm"] = np.log2(long_df["tpm"] + 1)
    long_df["z"] = long_df.groupby(["tissue", "gene"])["log_tpm"].transform(
        lambda x: (x - x.mean()) / x.std()
    )
    return long_df

mtdna_genes = ["MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
               "MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ATP6","MT-ATP8"]

oxphos_long = zscore_long(pd.read_csv("data/raw/gtex/gtex_oxphos_long.csv"))
oxphos_long = oxphos_long[~oxphos_long["gene"].isin(mtdna_genes)]  # nuclear only, as candidate predictors

candidates_long = zscore_long(pd.read_csv("data/raw/gtex/gtex_candidates_long.csv"))

features_long = pd.concat([oxphos_long[["sample_id","tissue","gene","z"]],
                            candidates_long[["sample_id","tissue","gene","z"]]])
features_wide = features_long.pivot_table(index=["sample_id","tissue"], columns="gene", values="z").reset_index()

target = pd.read_csv("results/phase2_signcorrected.csv")[["sample_id","tissue","regen_score_v3"]]

df = features_wide.merge(target, on=["sample_id","tissue"])
gene_cols = [c for c in features_wide.columns if c not in ("sample_id","tissue")]

candidate_genes = {"NAMPT", "SIRT1", "PPARGC1A"}

for tissue, sub in df.groupby("tissue"):
    sub = sub.dropna(subset=gene_cols + ["regen_score_v3"])
    if len(sub) < 50:
        continue
    X = sub[gene_cols]
    y = sub["regen_score_v3"]

    model = RandomForestRegressor(n_estimators=300, max_depth=6, random_state=42, n_jobs=-1)
    model.fit(X, y)

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)
    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    ranking = pd.Series(mean_abs_shap, index=gene_cols).sort_values(ascending=False)

    print(f"\n=== {tissue} (n={len(sub)}, R2={model.score(X, y):.3f}) ===")
    print("Top 15 genes by SHAP importance:")
    for i, (gene, val) in enumerate(ranking.head(15).items(), 1):
        flag = "  <-- CANDIDATE" if gene in candidate_genes else ""
        print(f"  {i}. {gene}: {val:.4f}{flag}")

    for cand in candidate_genes:
        if cand in ranking.index:
            rank = list(ranking.index).index(cand) + 1
            print(f"  {cand} overall rank: {rank}/{len(ranking)}")