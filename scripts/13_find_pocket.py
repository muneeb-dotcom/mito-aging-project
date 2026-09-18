import glob, os, json
import numpy as np
from Bio.PDB import PDBParser

SRC_DIR = "data/raw/pdb"
IGNORE_RESIDUES = {
    "HOH", "WAT", "NA", "CL", "MG", "ZN", "CA", "K", "MN", "FE", "CO", "NI",
    "SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "PEG", "PG4", "BME", "MPD",
    "TRS", "IMD", "FMT", "ACY", "1PE", "P6G", "EPE", "MES", "CIT", "MPO", "DTD",
    "MSE", "SEP", "TPO", "PTR", "CSO", "CSD", "MLY", "KCX", "PCA", "HYP",
    "NAG", "BMA", "MAN", "FUC", "GAL", "BGC", "GLC", "SIA", "NDG",
    "CLR", "CHS", "OLA", "OLB", "OLC", "LMT", "LDA", "LMN", "PLM", "MYR",
}

os.makedirs("data/pdbqt", exist_ok=True)
parser = PDBParser(QUIET=True)
configs = {}

for pdb_file in glob.glob(f"{SRC_DIR}/*.pdb"):
    name = os.path.basename(pdb_file).replace(".pdb", "")
    structure = parser.get_structure(name, pdb_file)

    # Each ligand instance kept SEPARATE (chain + resname + resnum), not merged
    instances = {}
    for residue in structure.get_residues():
        if residue.id[0] == " ":
            continue
        resname = residue.resname.strip()
        if resname in IGNORE_RESIDUES:
            continue
        chain_id = residue.get_parent().id
        key = (resname, chain_id, residue.id[1])
        coords = [atom.coord for atom in residue.get_atoms()]
        instances[key] = coords

    if not instances:
        print(name, "-- no candidate ligand")
        continue

    # Pick the biggest single ligand resname group (by number of instances with most atoms each)
    resname_counts = {}
    for (resname, chain_id, resnum), coords in instances.items():
        resname_counts.setdefault(resname, []).append((chain_id, resnum, coords))

    best_resname = max(resname_counts, key=lambda r: max(len(c) for _, _, c in resname_counts[r]))
    all_copies = resname_counts[best_resname]
    print(f"{name}: found {len(all_copies)} copies of {best_resname} (chains: {[c for c,_,_ in all_copies]})")

    # Use only the FIRST copy (one binding site) rather than merging all copies
    chain_id, resnum, coords = all_copies[0]
    coords = np.array(coords)
    center = coords.mean(axis=0)
    extent = coords.max(axis=0) - coords.min(axis=0)
    size = np.maximum(extent + 10, [20, 20, 20])

    if size.max() > 40:
        print(f"  SUSPICIOUS: box still too big ({size.round(1)}) — needs manual check")
        continue

    configs[name] = {
        "ligand_used": best_resname, "chain": chain_id,
        "center_x": round(float(center[0]), 2), "center_y": round(float(center[1]), 2), "center_z": round(float(center[2]), 2),
        "size_x": round(float(size[0]), 2), "size_y": round(float(size[1]), 2), "size_z": round(float(size[2]), 2),
    }
    print(f"  -> using chain {chain_id} copy: center={center.round(1)}, size={size.round(1)}")

with open("data/pdbqt/box_configs.json", "w") as f:
    json.dump(configs, f, indent=2)