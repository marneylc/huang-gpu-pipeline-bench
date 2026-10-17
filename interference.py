#!/usr/bin/env python
"""What else arrives at the same mass as acetyltropine?

Stage 5 of the slide. Being able to ionize a molecule is not the same as being able to see it.
Two kinds of neighbour matter in a MALDI imaging experiment:

  the matrix cluster series, which follows from the matrix formula alone and so costs nothing
  to work out, and changes depending on which matrix you pick;

  the handful of ions that dominate positive mode imaging of tissue whatever the matrix.

Each neighbour is reported with the resolving power needed to separate it, because whether an
overlap matters is a property of the instrument, not of the pair.
"""
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

HERE = Path(__file__).resolve().parent
PROTON, NA, K, H2O = 1.007276, 22.989218, 38.963158, 18.010565
NEAR_DA = 1.0

# the thirteen matrices from the screen, as formulas; the cluster masses follow from these
MATRICES = {
    "CHCA": "N#CC(=Cc1ccc(O)cc1)C(=O)O", "2,5-DHB": "OC(=O)c1cc(O)ccc1O",
    "3,5-DHB": "OC(=O)c1cc(O)cc(O)c1", "sinapinic": "COc1cc(/C=C/C(=O)O)cc(OC)c1O",
    "ferulic": "COc1cc(/C=C/C(=O)O)ccc1O", "9-AA": "Nc1c2ccccc2nc2ccccc12",
    "norharmane": "c1ccc2c(c1)[nH]c1cnccc12", "1,5-DAN": "Nc1cccc2c(N)cccc12",
    "DHAP": "CC(=O)c1ccc(O)cc1O", "THAP": "CC(=O)c1c(O)cc(O)cc1O",
    "dithranol": "O=C1c2ccccc2Cc2c1cccc2O", "MBT": "S=c1[nH]c2ccccc2s1",
    "PNA": "Nc1ccc([N+](=O)[O-])cc1",
}
# ions present in almost any tissue section in positive mode, whatever matrix you use
TISSUE = {"phosphocholine head group": 184.0733, "cholesterol minus water": 369.3516}


def cluster_ions(nmax=3):
    for name, smi in MATRICES.items():
        mw = ExactMolWt(Chem.MolFromSmiles(smi))
        for n in range(1, nmax + 1):
            yield f"{name} [{n}M+H]+", n * mw + PROTON
            yield f"{name} [{n}M+Na]+", n * mw + NA
            yield f"{name} [{n}M+K]+", n * mw + K
            yield f"{name} [{n}M-H2O+H]+", n * mw - H2O + PROTON


def main():
    smi = next(l.split()[0] for l in (HERE / "ligands.smi").read_text().splitlines()
               if l.strip() and not l.startswith("#") and l.split()[1] == "acetyltropine")
    mz = ExactMolWt(Chem.MolFromSmiles(smi)) + PROTON
    print(f"acetyltropine [M+H]+  {mz:.4f}\n")

    found = []
    for label, ion_mz in cluster_ions():
        gap = ion_mz - mz
        if abs(gap) <= NEAR_DA:
            found.append((abs(gap), label, ion_mz, gap, "matrix cluster"))
    for label, ion_mz in TISSUE.items():
        gap = ion_mz - mz
        if abs(gap) <= NEAR_DA:
            found.append((abs(gap), label, ion_mz, gap, "tissue background"))
    found.sort()

    if not found:
        print(f"nothing within {NEAR_DA} m/z.")
    else:
        print(f"{'neighbour':34s} {'m/z':>10s} {'gap':>10s} {'kind':>18s} {'R needed':>10s}")
        for _, label, ion_mz, gap, kind in found:
            print(f"{label:34s} {ion_mz:10.4f} {gap*1000:+9.0f} mDa {kind:>18s} "
                  f"{abs(mz/gap):10,.0f}")

    nearest_matrix = min((abs(m - mz) for _, m in cluster_ions()), default=None)
    print(f"\nnearest matrix cluster ion of any of the {len(MATRICES)} matrices: "
          f"{nearest_matrix:.2f} m/z away")
    print("So the matrices are not the problem here. The problem is the other column: the\n"
          "phosphocholine head group is the most abundant positive ion in almost any tissue\n"
          "section, and no choice of matrix moves it.")
    (HERE / "interference.json").write_text(json.dumps(
        {"mz": round(mz, 4),
         "neighbours": [{"ion": l, "mz": round(m, 4), "gap_mDa": round(g*1000, 1),
                         "kind": k, "resolving_power": round(abs(mz/g))}
                        for _, l, m, g, k in found]}, indent=1))


if __name__ == "__main__":
    main()
