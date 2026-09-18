import requests

pdb_url = "https://files.rcsb.org/download/2GVJ.pdb"
r = requests.get(pdb_url)
r.raise_for_status()
with open("data/raw/pdb/NAMPT_2GVJ.pdb", "wb") as f:
    f.write(r.content)
print("Saved NAMPT_2GVJ.pdb,", len(r.content), "bytes")

ligands = {
    "nicotinamide_CID936": 936,
    "NR_CID439924": 439924,
    "NMN_CID14180": 14180,
}

for name, cid in ligands.items():
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    r = requests.get(url)
    r.raise_for_status()
    with open(f"data/raw/{name}.sdf", "wb") as f:
        f.write(r.content)
    print(f"Saved {name}.sdf, {len(r.content)} bytes")