import sys, os, multiprocessing as mp, numpy as np
def prep_receptor():
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    lig=[l for l in open("6A93.pdb") if l.startswith("HETATM") and l[17:20]=="8NU" and l[21]=="A"]
    xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lig])
    c=xyz.mean(0); span=xyz.max(0)-xyz.min(0)
    open("box.txt","w").write(" ".join(f"{v:.3f}" for v in c)+"\n")
    print("ligand atoms",len(lig),"center",c.round(2),"span",span.round(1))
    prot=[l for l in open("6A93.pdb") if l.startswith("ATOM") and l[21]=="A"]+["END\n"]
    open("recA.pdb","w").writelines(prot)
    f=PDBFixer("recA.pdb"); f.findMissingResidues(); f.missingResidues={}
    f.findNonstandardResidues(); f.replaceNonstandardResidues(); f.removeHeterogens(False)
    f.findMissingAtoms(); f.addMissingAtoms(); f.addMissingHydrogens(7.4)
    PDBFile.writeFile(f.topology,f.positions,open("recA_H.pdb","w"))
    os.system("obabel recA_H.pdb -xr -O receptor.pdbqt 2>/dev/null")  # needs the dock env on PATH
def lig(args):
    smi,name=args
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    try:
        m=Chem.AddHs(Chem.MolFromSmiles(smi))
        if AllChem.EmbedMolecule(m,randomSeed=1)!=0: return None
        AllChem.MMFFOptimizeMolecule(m)
        setups=MoleculePreparation().prepare(m)
        s,ok,err=PDBQTWriterLegacy.write_string(setups[0])
        if not ok: return None
        open(f"ligs/{name}.pdbqt","w").write(s); return name
    except Exception: return None
if __name__=="__main__":
    prep_receptor()
    os.makedirs("ligs",exist_ok=True)
    pairs=[l.split() for l in open("../ligands_2000.smi") if l.strip()]
    with mp.Pool(int(os.environ.get("SLURM_CPUS_PER_TASK",8))) as p:
        ok=[r for r in p.map(lig,pairs,chunksize=20) if r]
    print("ligands prepared",len(ok),"of",len(pairs))
