import sys, os, time, numpy as np, torch, multiprocessing as mp
try:  # needed on ARM, where torch's NVRTC 12.6 rejects Warp's -pch option; harmless on x86
    import warp; warp.config.use_precompiled_headers = False
except ImportError:
    pass
from rdkit import Chem
from rdkit.Chem import AllChem
def confs(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi))
    ids=AllChem.EmbedMultipleConfs(m,numConfs=10,randomSeed=3)
    if len(ids)<10: return None
    return np.array([a.GetAtomicNum() for a in m.GetAtoms()]), np.array([m.GetConformer(i).GetPositions() for i in ids]), Chem.GetFormalCharge(m)
if __name__=="__main__":
    mode=sys.argv[1]
    if mode=="embed":
        smis=[l.split()[0] for l in open("../ligands_2000.smi")][:1000]
        with mp.Pool(int(os.environ.get("SLURM_CPUS_PER_TASK",8))) as p: res=[r for r in p.map(confs,smis,chunksize=10) if r]
        np.save("confs.npy", np.array(res, dtype=object), allow_pickle=True); print("molecules",len(res),"conformers",10*len(res)); sys.exit()
    from aimnet.calculators import AIMNet2Calculator
    res=np.load("confs.npy",allow_pickle=True)
    part=int(sys.argv[2]); nparts=int(sys.argv[3]); res=res[part::nparts]
    calc=AIMNet2Calculator("aimnet2-wb97m-d3_0", device="cuda")
    def ev(z,X,q):
        d={"coord":torch.tensor(X,dtype=torch.float32),"numbers":torch.tensor(np.tile(z,(len(X),1))),"charge":torch.full((len(X),),float(q))}
        return calc(d, forces=True)
    ev(*res[0]); torch.cuda.synchronize()
    t=time.time(); n=0; E=[]
    for z,X,q in res:
        out=ev(z,X,q); E.append(out["energy"].detach().cpu().numpy()); n+=len(X)
    torch.cuda.synchronize(); dt=time.time()-t
    print(f"part{part}/{nparts} conformers={n} wall_s={dt:.1f} conf_per_s={n/dt:.0f} E0={E[0][0]:.4f}", flush=True)
