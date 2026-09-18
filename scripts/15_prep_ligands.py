import glob, os, subprocess

os.makedirs("data/pdbqt/ligands", exist_ok=True)

for sdf in glob.glob("data/raw/*.sdf"):
    name = os.path.basename(sdf).replace(".sdf", "")
    out = f"data/pdbqt/ligands/{name}.pdbqt"
    r = subprocess.run(["obabel", sdf, "-O", out, "-p", "7.4", "--partialcharge", "gasteiger"],
                        capture_output=True, text=True)
    print(name, "OK" if os.path.exists(out) else "FAILED")