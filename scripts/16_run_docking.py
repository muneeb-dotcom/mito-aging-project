import os, json, re, subprocess
import pandas as pd

with open("data/pdbqt/box_configs.json") as f:
    box_configs = json.load(f)

PAIRS = [
    ("nicotinamide_CID936", "NAMPT_2GVJ"),
    ("NR_CID439924", "NAMPT_2GVJ"),
    ("NMN_CID14180", "NAMPT_2GVJ"),
]

VINA_EXE = r"C:\Users\CC\serum-primed-healing\vina.exe"
EXHAUSTIVENESS = 8
N_POSES = 9

os.makedirs("results/poses", exist_ok=True)
results = []

for ligand_name, target_name in PAIRS:
    receptor = f"data/pdbqt/receptors/{target_name}.pdbqt"
    ligand = f"data/pdbqt/ligands/{ligand_name}.pdbqt"
    box = box_configs.get(target_name)

    if not os.path.exists(receptor) or not os.path.exists(ligand) or box is None:
        print("SKIP", ligand_name, target_name)
        continue

    out_path = f"results/poses/{ligand_name}_{target_name}.pdbqt"
    cmd = [
        VINA_EXE,
        "--receptor", receptor, "--ligand", ligand,
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]), "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(N_POSES),
        "--out", out_path,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)

    score = None
    for line in r.stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.?\d*)", line)
        if m:
            score = float(m.group(1))
            break

    print(ligand_name, target_name, score)
    if score is None:
        print(r.stdout[-500:], r.stderr[-500:])
    results.append({"compound": ligand_name, "target": target_name, "docking_score_kcal_mol": score})

df = pd.DataFrame(results)
df.to_csv("results/phase4_docking_scores.csv", index=False)
print(f"\nSaved {len(df)} scores -> results/phase4_docking_scores.csv")
print(df)