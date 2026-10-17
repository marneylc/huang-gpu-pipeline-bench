#!/usr/bin/env python
"""Does the docking setup reproduce a pose we already know the answer to?

Tiotropium is the ligand in the 5CXV crystal structure, so its real position is known. Dock it
back in and measure how far the top pose lands from that position. Under 2.0 angstrom is the
usual pass mark. If a setup cannot reproduce a measured pose, its scores for everything else
mean nothing, and the scores will look perfectly reasonable while it fails.

RMSD is symmetry aware and computed in place, with no realignment: CalcRMS, not GetBestRMS.
Realigning would quietly forgive a pose that sits in the wrong place.
"""
import json, sys
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate
RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
CONTROL = "tiotropium"


def smiles_for(name):
    for line in (HERE / "ligands.smi").read_text().splitlines():
        if line.strip() and not line.startswith("#") and line.split()[1] == name:
            return line.split()[0]
    sys.exit(f"{name} not in ligands.smi")


def protonate(m):
    for a in m.GetAtoms():
        if (a.GetSymbol() == "N" and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 0
                and a.GetDegree() == 3 and not a.GetIsAromatic()):
            a.SetFormalCharge(1); a.SetNumExplicitHs(1); a.SetNoImplicit(True)
    Chem.SanitizeMol(m)
    return m


def main():
    smi = smiles_for(CONTROL)
    tmpl = protonate(Chem.MolFromSmiles(smi))
    pose_file = HERE / "poses" / f"{CONTROL}.pdbqt"
    if not pose_file.exists():
        sys.exit(f"{pose_file} not found; run dock.py first")

    ref = Chem.MolFromPDBFile(str(HERE / "xtal_ligand.pdb"), sanitize=False, removeHs=True,
                              proximityBonding=True)
    ref = AllChem.AssignBondOrdersFromTemplate(tmpl, ref)
    poses = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(str(pose_file), skip_typing=True))[0])
    rms = [round(rdMolAlign.CalcRMS(Chem.Mol(poses, confId=c.GetId()), ref), 2)
           for c in poses.GetConformers()]

    verdict = "PASS" if rms[0] < 2.0 else "FAIL"
    print(f"redock {CONTROL} into its own crystal structure (5CXV)")
    print(f"  top pose      {rms[0]:5.2f} A   {verdict} (pass mark 2.0)")
    print(f"  best of {len(rms)}     {min(rms):5.2f} A")
    (HERE / "redock.json").write_text(json.dumps(
        {"ligand": CONTROL, "rmsd_top": rms[0], "rmsd_best": min(rms), "rmsd_all": rms}, indent=1))


if __name__ == "__main__":
    main()
