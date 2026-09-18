import pandas as pd

url = "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls"

xls = pd.ExcelFile(url, engine="xlrd")
print("Sheet names:", xls.sheet_names)

for name in xls.sheet_names:
    df = pd.read_excel(xls, sheet_name=name, nrows=2)
    print(f"\n--- {name} ---")
    print(list(df.columns))