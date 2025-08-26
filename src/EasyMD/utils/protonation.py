import os
import subprocess
from pathlib import Path

# --- configuration ---
pdb_in = "input.pdb"
pdb_out_variants = "input_variants.pdb"
target_pH = 7.0
decision_buffer = 0.3  # only flip state if |pKa - pH| > buffer

# --- 1) run PROPKA to estimate pKa ---
def run_propka(pdb_path: str):
    """
    Returns a dict like:
      { ('ASP', 45, 'A'): 3.6, ('GLU', 12, 'A'): 4.5, ('HIS', 88, 'A'): 6.8, ... }
    Keys are (resname, resid, chainID).
    """
    try:
        import propka  # type: ignore
        # PROPKA's Python API differs across versions; CLI is the most robust path.
        raise ImportError  # force fallback to CLI for consistency below
    except Exception:
        # CLI fallback: creates <pdb>.pka file
        subprocess.run(["propka3", pdb_path], check=True)
        pka_file = Path(pdb_path).with_suffix(".pka")
        return parse_propka_pka(str(pka_file))

def parse_propka_pka(pka_file: str):
    pka_map = {}
    with open(pka_file) as f:
        for line in f:
            # PROPKA lines often look like:
            # " ASP A  45    pKa= 3.60 ..."
            # " HIS B  88    pKa= 6.80 ..."
            parts = line.strip().split()
            if len(parts) < 5: 
                continue
            res = parts[0]
            chain = parts[1]
            try:
                resid = int(parts[2])
            except ValueError:
                continue
            if "pKa=" in line:
                try:
                    pka_val = float(line.split("pKa=")[1].split()[0])
                    pka_map[(res, resid, chain)] = pka_val
                except Exception:
                    pass
    return pka_map

pka_map = run_propka(pdb_in)

# --- 2) decide protonation states from pKa vs pH ---
def decide_state(resname, pka, pH, buf=0.3):
    """
    Returns a variant tag or None if default state is fine.
    Logic:
      - ASP/GLU: protonated (ASH/GLH) if pKa > pH + buf
      - HIS: HIP if pKa > pH + buf; otherwise neutral tautomer (choose HID/HIE later)
      - CYS/TYR: protonated if pKa > pH + buf (rare at neutral pH unless special env.)
      - LYS/ARG: remain protonated unless pKa << pH - buf (rare)
    """
    acidic = {"ASP": ("ASH", 4.0), "GLU": ("GLH", 4.4), "TYR": (None, 10.1), "CYS": (None, 8.3)}
    basic  = {"LYS": (None, 10.4), "ARG": (None, 12.5), "HIS": (None, 6.5)}

    if resname in acidic:
        protonated_name, _ = acidic[resname]
        if pka is not None and pka > pH + buf:
            return protonated_name  # e.g., ASH/GLH; TYR/CYS have no Amber protonated alt name commonly used
        return None

    if resname == "HIS":
        if pka is not None and pka > pH + buf:
            return "HIP"  # double protonated (charged)
        # neutral tautomer choice (HID vs HIE) will be decided by H-bonding; default to HIE here
        return "HIE"

    if resname in basic:
        # Rarely deprotonated at neutral pH—keep default (protonated) unless very low pKa
        if pka is not None and pka < pH - buf:
            # If you truly need deprotonated LYS/ARG, you’d need specific residue names supported by your FF.
            return None
        return None

    return None

# Collect decisions
decisions = {}  # (chain,resid) -> new_resname
for (res, resid, chain), pka in pka_map.items():
    newname = decide_state(res, pka, target_pH, decision_buffer)
    if newname:
        decisions[(chain, resid)] = newname

# --- 3) rewrite residue names in a copy of the PDB ---
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure("prot", pdb_in)

for model in structure:
    for chain in model:
        for res in chain:
            het, resseq, icode = res.id
            if het != " ":  # skip HETATM/ligands here
                continue
            key = (chain.id, resseq)
            old = res.resname.strip()
            if (chain.id, resseq) in decisions:
                res.resname = decisions[key]
            elif old == "HIS" and (chain.id, resseq) not in decisions:
                # If HIS wasn’t in pKa map (e.g., old PROPKA output), choose a neutral tautomer
                res.resname = "HIE"

io = PDBIO()
io.set_structure(structure)
io.save(pdb_out_variants)

print(f"Saved residue-renamed PDB: {pdb_out_variants}")
"""
# --- 4) add hydrogens with OpenMM (optional but recommended) ---
# This step ensures correct H placement for the variants you just set.
from pdbfixer import PDBFixer
from openmm.app import Modeller, PDBFile, ForceField
from openmm import unit

fixer = PDBFixer(filename=pdb_out_variants)
fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms()

mod = Modeller(fixer.topology, fixer.positions)
ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")  # example
mod.addHydrogens(forcefield=ff, pH=target_pH)  # respects HID/HIE/HIP, ASH/GLH, etc.

with open("input_prepped.pdb", "w") as f:
    PDBFile.writeFile(mod.topology, mod.positions, f)
print("Wrote hydrogenated structure: input_prepped.pdb")
"""