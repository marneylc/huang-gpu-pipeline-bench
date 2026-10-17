#!/usr/bin/env python
"""Build the 1,997-compound screening library for the GPU scaling run.

Same preparation as prep.py, just a lot more of it, and in parallel because embedding and
minimising 2,000 conformers one at a time is the slow part.
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
OUT = HERE / "library"


def one(arg):
    smi, name = arg
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return name, "unparseable"
        # same protonation rule as prep.py, for the same reason
        for a in m.GetAtoms():
            if (a.GetSymbol() == "N" and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 0
                    and a.GetDegree() == 3 and not a.GetIsAromatic()):
                a.SetFormalCharge(1); a.SetNumExplicitHs(1); a.SetNoImplicit(True)
        Chem.SanitizeMol(m)
        mh = Chem.AddHs(m)
        if AllChem.EmbedMolecule(mh, randomSeed=1) != 0:
            return name, "embed failed"
        AllChem.MMFFOptimizeMolecule(mh)
        s, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mh)[0])
        if not ok:
            return name, str(err)[:40]
        (OUT / f"{name}.pdbqt").write_text(s)
        return name, None
    except Exception as e:
        return name, str(e)[:40]


def main():
    OUT.mkdir(exist_ok=True)
    jobs = []
    for line in (HERE / "ligands_2000.smi").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            p = line.split()
            jobs.append((p[0], p[1]))
    # the named tropanes go in too, so they can be ranked against the library
    for line in (HERE / "ligands.smi").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            p = line.split()
            jobs.append((p[0], p[1]))

    n = int(sys.argv[1]) if len(sys.argv) > 1 else None
    with ProcessPoolExecutor(n) as ex:
        res = list(ex.map(one, jobs, chunksize=16))
    bad = [(k, v) for k, v in res if v]
    print(f"prepared {len(res) - len(bad)} of {len(res)} ligands into {OUT.name}/")
    for k, v in bad[:10]:
        print(f"  skipped {k}: {v}")


if __name__ == "__main__":
    main()
