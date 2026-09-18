import pandas as pd

url = "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls"
df = pd.read_excel(url, sheet_name="A Human MitoCarta3.0", engine="xlrd")

symbol_col = "Symbol"
pathway_col = "MitoCarta3.0_MitoPathways"

oxphos_mask = df[pathway_col].astype(str).str.contains("OXPHOS", case=False, na=False)
nuclear_oxphos = df.loc[oxphos_mask, symbol_col].dropna().unique().tolist()

# Drop the 13 mtDNA-encoded genes (MT-*) — tracked separately as the mtDNA gene set
nuclear_oxphos = [g for g in nuclear_oxphos if not g.startswith("MT-")]

print(f"Nuclear-encoded OXPHOS genes found: {len(nuclear_oxphos)}")
print(nuclear_oxphos[:10], "...")

out = pd.DataFrame({"gene_symbol": nuclear_oxphos})
out.to_csv("data/raw/mitocarta/nuclear_oxphos_genes.csv", index=False)
print("Saved to data/raw/mitocarta/nuclear_oxphos_genes.csv")