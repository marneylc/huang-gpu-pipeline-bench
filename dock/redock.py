import sys, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
SMI="Cc1nc2n(c(=O)c1CCN1CCC(CC1)c1noc3cc(F)ccc13)CCCC2"
if sys.argv[1]=="prep":
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    m=Chem.AddHs(Chem.MolFromSmiles("Cc1nc2n(c(=O)c1CC[NH+]1CCC(CC1)c1noc3cc(F)ccc13)CCCC2"))
    AllChem.EmbedMolecule(m,randomSeed=1); AllChem.MMFFOptimizeMolecule(m)
    s,ok,_=PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0]); open("risperidone.pdbqt","w").write(s)
    lines=[l for l in open("6A93.pdb") if l.startswith("HETATM") and l[17:20]=="8NU" and l[21]=="A"]
    open("xtal_lig.pdb","w").writelines(lines+["END\n"])
else:
    from meeko import PDBQTMolecule, RDKitMolCreate
    ref=Chem.MolFromPDBFile("xtal_lig.pdb"); ref=AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(SMI),ref)
    pm=PDBQTMolecule.from_file("redock/risperidone_out.pdbqt", skip_typing=True)
    pose=RDKitMolCreate.from_pdbqt_mol(pm)[0]
    pose=Chem.RemoveHs(pose); pose=Chem.MolFromSmiles(Chem.MolToSmiles(pose)) and pose
    score=[l for l in open("redock/risperidone_out.pdbqt") if "VINA RESULT" in l][0].split()[3]
    rms=rdMolAlign.CalcRMS(Chem.RemoveHs(pose,sanitize=False),ref)
    print(f"risperidone redock: top score {score} kcal/mol, heavy-atom RMSD to crystal {rms:.2f} A")
