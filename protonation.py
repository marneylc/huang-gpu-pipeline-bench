#!/usr/bin/env python
"""What does getting the protonation wrong actually cost?

Every tropane here is a tertiary amine with a pKa near 9.5, so at pH 7.4 it is more than 99
percent protonated, and that ammonium makes the salt bridge that anchors it in the pocket.
prep.py builds them that way. This script builds them the other way as well, changes nothing
else, docks both, and puts the two side by side.

Read the two right-hand columns together. The score barely moves. The pose sometimes moves
several angstrom. Nothing in the score tells you which case you are in, and a screen that
reports only scores cannot distinguish them.
"""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from vina import Vina
RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
BOX, EXHAUSTIVENESS, N_POSES = 22.0, 32, 9


def build_mol(smi, protonated):
    """The molecule with its tertiary amines protonated, or not."""
    m = Chem.MolFromSmiles(smi)
    if protonated:
        for a in m.GetAtoms():
            if (a.GetSymbol() == "N" and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 0
                    and a.GetDegree() == 3 and not a.GetIsAromatic()):
                a.SetFormalCharge(1); a.SetNumExplicitHs(1); a.SetNoImplicit(True)
        Chem.SanitizeMol(m)
    return m


def build(smi, protonated):
    mh = Chem.AddHs(build_mol(smi, protonated))
    AllChem.EmbedMolecule(mh, randomSeed=1); AllChem.MMFFOptimizeMolecule(mh)
    s, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mh)[0])
    assert ok, err
    return s


def top_pose(path):
    """Top pose, charges stripped so the two charge states are comparable as heavy-atom frames."""
    m = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(path), skip_typing=True))[0]
    m = Chem.RemoveHs(Chem.Mol(m, confId=m.GetConformers()[0].GetId()))
    for a in m.GetAtoms():
        a.SetFormalCharge(0); a.SetNumExplicitHs(0); a.SetNoImplicit(False)
    Chem.SanitizeMol(m, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    return m


def main():
    cx, cy, cz = [float(x) for x in (HERE / "box.txt").read_text().split()]
    v = Vina(sf_name="vina", seed=42, verbosity=0)
    v.set_receptor(str(HERE / "receptor.pdbqt"))
    v.compute_vina_maps(center=[cx, cy, cz], box_size=[BOX] * 3)

    work = HERE / "protonation"; work.mkdir(exist_ok=True)
    rows = {}
    for line in (HERE / "ligands.smi").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        smi, name = line.split()[:2]
        got = {}
        for state in ("protonated", "neutral"):
            f = work / f"{name}_{state}.pdbqt"
            f.write_text(build(smi, state == "protonated"))
            v.set_ligand_from_file(str(f))
            v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
            out = work / f"{name}_{state}_out.pdbqt"
            v.write_poses(str(out), n_poses=N_POSES, overwrite=True)
            got[state] = (float(np.array(v.energies(n_poses=N_POSES))[0, 0]), out)
        # a quaternary ammonium like tiotropium is already charged, so protonate() cannot
        # change it and the two builds are the same molecule. Say so rather than printing a
        # zero that looks like a result.
        same = Chem.MolToSmiles(build_mol(smi, False)) == Chem.MolToSmiles(build_mol(smi, True))
        try:
            shift = None if same else round(rdMolAlign.CalcRMS(
                top_pose(got["protonated"][1]), top_pose(got["neutral"][1])), 2)
        except Exception:
            shift = None
        rows[name] = {"score_protonated": round(got["protonated"][0], 2),
                      "score_neutral": round(got["neutral"][0], 2),
                      "score_change": round(got["protonated"][0] - got["neutral"][0], 2),
                      "pose_shift_A": shift}

    (HERE / "protonation.json").write_text(json.dumps(rows, indent=1))
    print(f"{'ligand':24s} {'protonated':>11s} {'neutral':>9s} {'score moves':>12s} {'pose moves':>11s}")
    for n, r in sorted(rows.items(), key=lambda t: -(t[1]["pose_shift_A"] or 0)):
        sh = "no basic N" if r["pose_shift_A"] is None else f"{r['pose_shift_A']:.2f} A"
        print(f"{n:24s} {r['score_protonated']:11.2f} {r['score_neutral']:9.2f} "
              f"{r['score_change']:+12.2f} {sh:>11s}")
    moved = [r for r in rows.values() if r["pose_shift_A"]]
    if moved:
        print(f"\nlargest score change {max(abs(r['score_change']) for r in moved):.2f} kcal/mol, "
              f"largest pose shift {max(r['pose_shift_A'] for r in moved):.2f} A")


if __name__ == "__main__":
    main()
