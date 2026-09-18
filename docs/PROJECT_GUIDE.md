\# Mitonuclear Imbalance \& Aging — Full Project Guide



A walkthrough of what we asked, what we built, what the data actually said,

and what's next — written so both a non-biologist and a biologist can follow

every step.



\## Project overview



\*\*In plain terms:\*\* Every cell has tiny power plants called mitochondria.

Building those power plants needs instructions from two different places:

most instructions live in the cell's main DNA library (the nucleus), but 13

instructions live inside the mitochondria itself, in their own small separate

DNA. As we age, these two instruction sets can drift out of sync. This

project asked: does that drift happen with age? Does it make tissue worse at

healing/regrowing? Is there one gene we could target with a drug to fix it?



\*\*In scientific terms:\*\* Oxidative phosphorylation (OXPHOS) complexes are

built from proteins encoded by two separate genomes: \~1,100+ nuclear genes

and 13 mitochondrial-DNA (mtDNA) genes, which must be expressed in tight

stoichiometric balance. \*\*Mitonuclear imbalance\*\* — a drift in the ratio of

nuclear-to-mtDNA OXPHOS gene expression — is a documented feature of aging.

This project tested whether that imbalance (1) increases with age across

human tissues, (2) predicts suppression of a pre-established regeneration

gene signature, (3) is driven by an identifiable chokepoint gene, and (4) is

druggable — then attempted (5) a cross-species validation in a short-lived

fish model.



\---



\## Phase 0 — Environment setup ✅



\*\*In plain terms:\*\* Like setting up a kitchen before cooking — a

self-contained software workspace with pre-built biology packages installed.



\*\*In scientific terms:\*\* Isolated Conda environment (`mito-aging`, Python

3.10 + R 4.3, supplemented by the machine's R 4.6.1 for Bioconductor) with

the full Python/R stack for RNA-seq retrieval, statistics, ML, and docking.



\*\*Tools:\*\*

| Tool | What it is | Why used | Manual alternative |

|---|---|---|---|

| Conda | Package/environment manager | Keeps project dependencies isolated | Install globally (risky, version clashes) |

| Python + pip | Programming language + package installer | Runs analysis, ML, docking scripts | No practical alternative |

| R + Bioconductor | Statistics language + biology package library | Some genomics tools (recount3, DESeq2) only exist here | No manual equivalent |



\*\*Outcome:\*\* working environment, one limitation — AutoDock Vina has no

native Windows Conda build, worked around by reusing a standalone `vina.exe`

from an earlier project.



\---



\## Phase 1 — Building the mitonuclear score ✅



\*\*In plain terms:\*\* Built a list of genes from each "instruction set" (13

mitochondrial, always the same; \~156 nuclear ones from a public catalog),

measured how active each gene was across real donated tissue samples, and

divided nuclear activity by mitochondrial activity to get one number per

sample. Higher = more imbalance.



\*\*In scientific terms:\*\* Filtered MitoCarta 3.0 for OXPHOS-pathway genes,

excluding the 13 mtDNA-encoded genes → 156 nuclear OXPHOS genes. Pulled bulk

RNA-seq (TPM, via recount3) for 5 GTEx tissues, computed

`mitonuclear\_score = mean(nuclear OXPHOS TPM) / mean(mtDNA OXPHOS TPM)` per

sample, tested Spearman correlation with donor age.



\*\*Tools:\*\*

| Tool | What it does | Outcome | Manual alternative |

|---|---|---|---|

| MitoCarta 3.0 | Broad Institute's curated mitochondrial gene catalog | 156 genes identified, 132 found expressed in GTEx | Manually filter the spreadsheet by eye |

| pandas + xlrd | Reads legacy `.xls` files | Clean CSV gene list exported | Filter in Excel |

| recount3 (R) | Public archive of reprocessed RNA-seq incl. GTEx | 1,083,420 rows across 5 tissues (6,945 samples) | Download raw FASTQ and align yourself |

| recount::getTPM() | Converts coverage to TPM (normalized expression) | Comparable expression across samples | Normalize manually by gene length/library size |

| SciPy Spearman correlation | Tests if two variables rise/fall together | See results below | Rank and compute by hand |



\*\*Results:\*\* 132/156 nuclear genes + 13/13 mtDNA genes found expressed (145

unique genes total).



| Tissue | Spearman ρ (score vs. age) | p-value | n |

|---|---:|---:|---:|

| BRAIN | −0.141 | 1.4e-14 | 2,931 |

| HEART | −0.094 | 0.0040 | 942 |

| LIVER | −0.164 | 0.0093 | 251 |

| MUSCLE | +0.184 | 3.8e-8 | 881 |

| SKIN | +0.084 | 0.00023 | 1,940 |



\*\*Outcome:\*\* tissue-specific, not universal — rises with age in muscle/skin,

falls in brain/heart/liver. All significant but modest effect sizes (|ρ|<0.2).



\---



\## Phase 2 — Testing the regeneration link ✅



\*\*In plain terms:\*\* Using a gene list from an earlier project that marks

active healing vs. scarring, checked whether high mitonuclear imbalance lines

up with \*less\* activity from the "good healing" genes. This phase is also a

case study in catching two real analysis mistakes before trusting a result.



\*\*In scientific terms:\*\* Pulled GTEx expression for the 111-gene

regeneration signature (4-gene core + 108-gene ECM/PI3K-Akt set) from a prior

DESeq2 mouse regeneration-vs-scarring analysis, correlated against the

mitonuclear score.



\*\*Tools:\*\*

| Tool | What it does | Why used |

|---|---|---|

| Prior DESeq2 output | Differential expression (log2FoldChange) from an earlier project | Gives the true direction each gene should move if pro-regenerative |

| NumPy log + z-score normalization | Rescales each gene to a comparable 0-centered scale | Stops high-expression genes dominating the average |

| statsmodels (OLS residuals, mixed-effects) | Regression / confound control / pooling across groups | Rules out RNA-quality confounds, gives one pooled estimate |



\*\*Mistake 1 — scale artifact:\*\* raw correlation was strongly positive

everywhere (ρ 0.44–0.84) — suspiciously uniform. After proper normalization,

SKIN's correlation collapsed to \~0, confirming part of the effect was just

unnormalized expression scale.



\*\*Mistake 2 — no sign correction:\*\* genes were averaged as if "more

expression = more regeneration" for all of them, but some are scarring

markers. Using the real DESeq2 fold-change direction to flip those genes

reversed the correlation in most tissues — matching what the hypothesis

actually predicted.



\*\*Final results\*\* (partial ρ, controlling for RNA quality/ischemic time):



| Tissue | Partial ρ | Direction |

|---|---:|---|

| BRAIN | −0.433 | Supports hypothesis |

| LIVER | −0.184 | Supports hypothesis |

| HEART | −0.124 | Weakly supports |

| MUSCLE | −0.004 | No effect |

| SKIN | +0.141 | Contradicts hypothesis |



\*\*Pooled mixed-effects model\*\* (all 6,926 samples, tissue as random effect):

`mitonuclear\_score` coefficient = −0.015, p<0.001 — supports the hypothesis

on average, though heavily weighted by BRAIN (42% of pooled n); the

tissue-by-tissue table is the more honest summary.



\*\*Why SKIN likely disagrees:\*\* the regeneration signature came from

comparing actively injured/regenerating tissue to scarring tissue in mice.

GTEx skin samples are normal, uninjured, resting-state tissue — an

injury-response signature may simply not transfer there.



\---



\## Phase 3 — Finding the chokepoint gene ✅



\*\*In plain terms:\*\* Trained a model to predict "how regenerative" a sample

looks just from 135 genes' activity, then asked which genes it actually

relied on most (via SHAP — like a judge explaining which evidence mattered).



\*\*In scientific terms:\*\* Per tissue, trained a Random Forest regressor (135

features: 132 nuclear OXPHOS genes + NAMPT/SIRT1/PPARGC1A) to predict the

sign-corrected regeneration score, ranked features by SHAP importance.



\*\*Tools:\*\*

| Tool | What it does | Manual alternative |

|---|---|---|

| scikit-learn RandomForestRegressor | ML model averaging many decision trees | Fit linear regression by hand (misses non-linear effects) |

| SHAP | Splits a prediction's "credit" fairly across input features | Drop one gene at a time and measure impact (slow, less rigorous) |



\*\*Results:\*\*



| Tissue | R² | Top predictor | NAMPT rank/135 | SIRT1 rank | PPARGC1A rank |

|---|---:|---|---:|---:|---:|

| BRAIN | 0.765 | COX20 | 9 | 43 | 39 |

| HEART | 0.741 | NDUFS3 | 2 | 29 | 11 |

| LIVER | 0.834 | PNKD | 7 | 48 | 26 |

| MUSCLE | 0.584 | CMC2 | 4 | 23 | 6 |

| SKIN | 0.885 | COX4I2 | 14 | 124 | 2 |



\*\*Outcome:\*\* "NAMPT is THE top predictor" wasn't confirmed — the top gene

differs by tissue, and each is a mitochondrial complex \*assembly factor\*

(helper protein, not a structural subunit). NAMPT's real finding is

consistency: the only candidate ranking top-15 in all five tissues, even

without ever being #1.



\---



\## Phase 4 — Docking candidate drugs ✅



\*\*In plain terms:\*\* Proteins have pockets shaped like a lock for a specific

key. Docking software tests many ways a molecule could sit in that pocket

and scores the fit — more negative = stronger fit. Tested three NAD⁺-related

supplements against NAMPT's pocket.



\*\*In scientific terms:\*\* Docked nicotinamide, NMN, and NR against the NAMPT

active site (PDB 2GVJ, co-crystallized with FK866, which competes directly

with the natural nicotinamide substrate) using AutoDock Vina, pocket defined

from FK866's crystallographic coordinates.



\*\*Tools:\*\*

| Tool | What it does | Manual alternative |

|---|---|---|

| RCSB PDB | Public database of solved 3D protein structures | Solve the structure yourself via X-ray crystallography |

| PubChem | Public database of chemical compound 3D structures | Build 3D coordinates by hand |

| Biopython (PDBParser) | Reads/writes protein structure files | Open in a viewer and read off coordinates manually |

| Open Babel | Converts chemistry file formats, assigns partial charges | No practical manual equivalent |

| AutoDock Vina | Docking engine — scores many poses by estimated binding energy | Visual inspection only (no numeric score) |



\*\*Results:\*\*



| Compound | Docking score (kcal/mol) | Biological role at NAMPT |

|---|---:|---|

| Nicotinamide | −6.704 | NAMPT's actual natural substrate |

| NMN | −7.528 | NAMPT's reaction product |

| NR | −8.459 | Not a NAMPT substrate — made by a different enzyme |



\*\*Caveat that matters:\*\* nicotinamide is literally what NAMPT converts into

NMN, so it should dock well — this validates the pipeline. NR is converted to

NMN by a different enzyme (nicotinamide riboside kinase), not NAMPT — its top

score likely reflects structural resemblance rather than real engagement.



\---



\## Phase 5 — Cross-species check ⚠️ Attempted, blocked



\*\*In plain terms:\*\* Tried to repeat Phase 1 in a fast-aging fish (the

turquoise killifish). Found two good public datasets, but both were missing

mitochondrial gene data entirely — a structural issue, not bad luck.



\*\*In scientific terms:\*\* Evaluated two independent \*Nothobranchius furzeri\*

RNA-seq resources: GSE66712 (JenAge, 5 tissues, numeric ages) and GSE308970

(2026 Nature Aging atlas, 13 tissues, 677 samples). Both had solid nuclear

OXPHOS ortholog coverage but zero mtDNA gene annotations in their processed

matrices.



\*\*Tools:\*\*

| Tool | What it does |

|---|---|

| GEOquery (R/Bioconductor) | Pulls public expression datasets/metadata from NCBI GEO |

| curl (resumable download, `-C -`) | Resumed a large (\~120MB) download after repeated connection drops instead of restarting from zero |



\*\*Why it's blocked:\*\* the standard RNA-seq quantification pipeline for

killifish (align → count against a gene annotation file) doesn't include the

mitochondrial genome as an annotated feature in either dataset checked — the

mtDNA reference is likely deposited/annotated separately from the nuclear

genome assembly these pipelines use. A structural gap in how killifish

transcriptomes are typically published, not a fixable data-selection mistake.



\---



\## Full results summary



| Phase | Question | Headline result |

|---|---|---|

| 1 | Does imbalance rise with age? | Tissue-specific: up in muscle/skin, down in brain/heart/liver |

| 2 | Does imbalance suppress regeneration genes? | Yes in brain/liver/heart (ρ up to −0.43), no in muscle, reversed in skin |

| 3 | Is there one chokepoint gene? | No single universal gene, but NAMPT is consistently top-15 in all 5 tissues |

| 4 | Is NAMPT druggable? | Nicotinamide and NMN dock favorably and are mechanistically valid; NR's top score is likely an artifact |

| 5 | Does the pattern hold in a fast-aging fish? | Untestable with current public data |



\*\*Defensible headline claim:\*\* age-related mitonuclear OXPHOS imbalance is

associated with reduced regeneration-signature expression in specific

tissues (most strongly brain and liver), driven less by a single master gene

than by tissue-specific OXPHOS assembly-factor bottlenecks, with NAMPT

standing out as the most consistent cross-tissue signal and a plausible,

druggable point of intervention.



\---



\## Future directions



\### Quick wins (days, same data already in hand)

\- Re-run the Phase 3 SHAP analysis restricted to BRAIN and LIVER only (the

&#x20; two tissues where Phase 2 actually held up) for a cleaner chokepoint ranking.

\- Dock a wider panel of known NAD⁺ pathway modulators (FK866 itself as a

&#x20; positive control, P7C3 compounds) to calibrate what a "strong" score means here.

\- Write this up as a short preprint/lab report — the data and figures already

&#x20; support a results section.



\### Moderate extensions (1–2 weeks, some new data)

\- Repeat Phase 2 with additional GTEx tissues to see how many fall on the

&#x20; "supports" side versus "contradicts."

\- Test whether the actual Phase 3 top chokepoint genes (COX20, NDUFS3, PNKD,

&#x20; CMC2, COX4I2) themselves change with age — only the 3 pre-registered

&#x20; candidates were tested for age-correlation, not the real top hits.

\- Re-attempt Phase 5 with zebrafish (Danio rerio), which has more mature,

&#x20; mitochondrially-complete genome annotation than killifish.



\### Full-scale extension (months, new experiments)

\- Wet-lab validation: measure NAMPT protein levels and NAD⁺ concentration

&#x20; directly in aged vs. young tissue (brain/liver especially).

\- Realign raw killifish (or axolotl/zebrafish) FASTQ data against a genome

&#x20; build retaining full mtDNA annotation, to properly complete Phase 5.

\- Test NAD⁺ precursor compounds experimentally in a regeneration model

&#x20; (zebrafish fin clip or mouse digit-tip regeneration) to see if

&#x20; supplementation improves regenerative outcomes.

\- Build a single-cell resolution version of this analysis — bulk tissue

&#x20; averages could be masking cell-type-composition shifts rather than true

&#x20; gene regulation changes.

