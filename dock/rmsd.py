from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate
SMI="Cc1nc2n(c(=O)c1CCN1CCC(CC1)c1noc3cc(F)ccc13)CCCC2"
tmpl=Chem.MolFromSmiles(SMI)
ref=AllChem.AssignBondOrdersFromTemplate(tmpl, Chem.MolFromPDBFile("xtal_lig.pdb"))
pm=PDBQTMolecule.from_file("redock/risperidone_out.pdbqt", skip_typing=True)
for i,pose in enumerate(RDKitMolCreate.from_pdbqt_mol(pm)[:1] and [None]):
    pass
mols=RDKitMolCreate.from_pdbqt_mol(pm)
pose=mols[0]
scores=[float(l.split()[3]) for l in open("redock/risperidone_out.pdbqt") if "VINA RESULT" in l]
best=None
for k in range(pose.GetNumConformers()):
    p=Chem.Mol(pose, confId=k) if False else Chem.Mol(pose)
    pc=Chem.RemoveHs(p)
    for a in pc.GetAtoms(): a.SetFormalCharge(0); a.SetNumExplicitHs(0)
    pc=AllChem.AssignBondOrdersFromTemplate(tmpl, Chem.RWMol(pc).GetMol())
    r=rdMolAlign.CalcRMS(pc, ref, prbId=k)
    print(f"pose {k+1} score {scores[k]:.2f} RMSD {r:.2f} A")
