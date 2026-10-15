import sys, time, warnings; warnings.filterwarnings("ignore")
from pyscf import gto
from gpu4pyscf.dft import rks, uks
xyz, chg, mult, tag = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
lines = open(xyz).read().splitlines()[2:]
m = gto.M(atom="\n".join(lines), basis="def2-tzvp", charge=chg, spin=mult-1, verbose=0)
mf = (rks.RKS(m, xc="HYB_GGA_XC_WB97X_D3") if mult == 1 else uks.UKS(m, xc="HYB_GGA_XC_WB97X_D3")).density_fit()
t = time.time(); e = mf.kernel()
print(f"{tag} {xyz} natm={m.natm} nao={m.nao} E={e:.8f} conv={mf.converged} wall_s={time.time()-t:.1f}", flush=True)
