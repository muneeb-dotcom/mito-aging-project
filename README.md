\# Mitonuclear Imbalance \& Aging



Testing whether the ratio of nuclear-to-mtDNA OXPHOS gene expression rises with

age across human tissues, whether that imbalance predicts suppression of a

regeneration gene signature, whether one "chokepoint" gene drives it, and

whether that gene is druggable.



\## Summary of results



| Phase | Question | Headline result |

|---|---|---|

| 0 | Environment setup | Done |

| 1 | Does imbalance rise with age? | Tissue-specific — up in muscle/skin, down in brain/heart/liver |

| 2 | Does imbalance suppress regeneration genes? | Yes in brain/liver/heart (partial ρ up to −0.43), no effect in muscle, reversed in skin |

| 3 | Is there one chokepoint gene? | No single universal gene; NAMPT is the only candidate ranking top-15/135 in all 5 tissues |

| 4 | Is NAMPT druggable with NAD⁺ precursors? | Nicotinamide (−6.70 kcal/mol) and NMN (−7.53) dock favorably and are mechanistically valid; NR (−8.46) is likely a structural-mimicry false positive |

| 5 | Cross-species check (killifish) | Blocked — both public datasets checked lack mtDNA gene quantification |



Full write-up with dual-tier (plain-language + technical) explanations, tool

rationale, and future directions: see \[`docs/PROJECT\_GUIDE.md`](docs/PROJECT\_GUIDE.md).



\## Repo structure

scripts/ all analysis scripts, in run order (01\_... through 23\_...)

results/ output tables (scores, docking results, merged datasets)

figures/ plots

data/ raw/intermediate data (not pushed — see .gitignore; scripts regenerate it)

docs/ full project guide


\## Reproducing



Environment: Conda (`mito-aging` env, Python 3.10 + R 4.3/4.6), see scripts for

package list. Scripts are numbered in intended run order; each `NN\_\*.py`/`.R`

file corresponds to one step described in `docs/PROJECT\_GUIDE.md`.



\## Tech stack



Python (pandas, scipy, scikit-learn, shap, statsmodels, rdkit, meeko),

R/Bioconductor (recount3, recount, DESeq2, biomaRt, tidyverse), AutoDock Vina,

Open Babel, MitoCarta 3.0, GTEx (via recount3), GEO.

