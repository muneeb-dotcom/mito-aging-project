import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv("results/phase2_signcorrected.csv").dropna(
    subset=["mitonuclear_score_v2", "regen_score_v3", "rin", "ischemic_time"]
)

model = smf.mixedlm(
    "regen_score_v3 ~ mitonuclear_score_v2 + rin + ischemic_time",
    data=df,
    groups=df["tissue"]
)
result = model.fit()
print(result.summary())