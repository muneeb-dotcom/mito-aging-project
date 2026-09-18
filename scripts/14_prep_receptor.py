import os
from Bio.PDB import PDBParser, PDBIO, Select

class ProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "

os.makedirs("data/pdb_clean", exist_ok=True)
os.makedirs("data/pdbqt/receptors", exist_ok=True)

parser = PDBParser(QUIET=True)
io = PDBIO()

name = "NAMPT_2GVJ"
structure = parser.get_structure(name, f"data/raw/pdb/{name}.pdb")
clean_path = f"data/pdb_clean/{name}_clean.pdb"
io.set_structure(structure)
io.save(clean_path, ProteinOnly())

import subprocess
out_path = f"data/pdbqt/receptors/{name}.pdbqt"
r = subprocess.run(["obabel", clean_path, "-O", out_path, "-xr", "-p", "7.4"],
                    capture_output=True, text=True)
print(name, "OK" if os.path.exists(out_path) else "FAILED")
print(r.stdout[-500:], r.stderr[-500:])