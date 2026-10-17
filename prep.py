#!/usr/bin/env python
"""Prepare the M1 muscarinic receptor and the tropane ligands for docking.

Receptor: PDB 5CXV, chain A, with tiotropium (ligand code 0HK) bound. Tiotropium is a tropane,
so redocking it checks the setup on the same scaffold as the compounds we care about.

The ligands are built with the ring nitrogen PROTONATED. Every tropane here is a tertiary amine
with a pKa near 9.5, so at pH 7.4 more than 99 percent of it carries a proton, and that ammonium
is what forms the salt bridge to Asp105 that holds the ligand in place. Building them neutral
removes the hydrogen bond donor the anchor needs. redock.py measures what that costs.
"""
import os, sys, urllib.request
from pathlib import Path
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

HERE = Path(__file__).resolve().parent
PDB_ID, CHAIN, XTAL_LIG = "5CXV", "A", "0HK"


def protonate_amines(m):
    """Put a proton on every non-aromatic tertiary amine, and report the resulting charge."""
    for a in m.GetAtoms():
        if (a.GetSymbol() == "N" and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 0
                and a.GetDegree() == 3 and not a.GetIsAromatic()):
            a.SetFormalCharge(1); a.SetNumExplicitHs(1); a.SetNoImplicit(True)
    Chem.SanitizeMol(m)
    return m, Chem.GetFormalCharge(m)


def main():
    pdb = HERE / f"{PDB_ID}.pdb"
    if not pdb.exists():
        print(f"downloading {PDB_ID}")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{PDB_ID}.pdb", pdb)

    # the crystal ligand gives both the docking box centre and the reference pose for redock.py
    L = [l for l in open(pdb) if l.startswith("HETATM") and l[17:20] == XTAL_LIG and l[21] == CHAIN]
    if not L:
        sys.exit(f"no {XTAL_LIG} found in chain {CHAIN} of {PDB_ID}")
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in L])
    (HERE / "box.txt").write_text(" ".join(f"{v:.3f}" for v in xyz.mean(0)) + "\n")
    (HERE / "xtal_ligand.pdb").write_text("".join(L))

    (HERE / "receptor_raw.pdb").write_text(
        "".join(l for l in open(pdb) if l.startswith("ATOM") and l[21] == CHAIN) + "END\n")
    f = PDBFixer(str(HERE / "receptor_raw.pdb"))
    f.findMissingResidues(); f.missingResidues = {}
    f.findNonstandardResidues(); f.replaceNonstandardResidues()
    f.removeHeterogens(False)
    f.findMissingAtoms(); f.addMissingAtoms(); f.addMissingHydrogens(7.4)
    PDBFile.writeFile(f.topology, f.positions, open(HERE / "receptor_H.pdb", "w"))
    os.system(f"obabel {HERE/'receptor_H.pdb'} -xr -O {HERE/'receptor.pdbqt'} 2>/dev/null")
    print(f"receptor {PDB_ID} chain {CHAIN}, box centre {xyz.mean(0).round(2)}")

    ligs = HERE / "ligands"; ligs.mkdir(exist_ok=True)
    for line in (HERE / "ligands.smi").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        smi, name = line.split()[:2]
        m, q = protonate_amines(Chem.MolFromSmiles(smi))
        mh = Chem.AddHs(m)
        AllChem.EmbedMolecule(mh, randomSeed=1)
        AllChem.MMFFOptimizeMolecule(mh)
        s, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mh)[0])
        if not ok:
            sys.exit(f"{name}: {err}")
        (ligs / f"{name}.pdbqt").write_text(s)
        print(f"  {name:22s} {m.GetNumHeavyAtoms():3d} heavy atoms  charge {q:+d}")


if __name__ == "__main__":
    main()
