#!/usr/bin/env python
"""Dock the tropane ligands into the M1 muscarinic receptor with AutoDock Vina.

Reports two columns, and the difference between them is the point. Raw score rewards size
almost automatically, because a bigger molecule makes more contacts. Score per heavy atom does
not. A small ligand can look poor on one and good on the other, so quoting only one of them
decides the answer before you have looked at it.
"""
import json, sys
from pathlib import Path
import numpy as np
from vina import Vina

HERE = Path(__file__).resolve().parent
EXHAUSTIVENESS, N_POSES, BOX = 32, 9, 22.0


def main():
    cx, cy, cz = [float(x) for x in (HERE / "box.txt").read_text().split()]
    cpu = int(sys.argv[1]) if len(sys.argv) > 1 else 0   # 0 lets Vina choose
    v = Vina(sf_name="vina", cpu=cpu, seed=42, verbosity=0)
    v.set_receptor(str(HERE / "receptor.pdbqt"))
    v.compute_vina_maps(center=[cx, cy, cz], box_size=[BOX] * 3)

    out_dir = HERE / "poses"; out_dir.mkdir(exist_ok=True)
    rows = {}
    for p in sorted((HERE / "ligands").glob("*.pdbqt")):
        heavy = sum(1 for l in p.read_text().splitlines()
                    if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("H", "HD"))
        v.set_ligand_from_file(str(p))
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        e = np.array(v.energies(n_poses=N_POSES))
        v.write_poses(str(out_dir / p.name), n_poses=N_POSES, overwrite=True)
        best = float(e[0, 0])
        rows[p.stem] = {"score": round(best, 3), "heavy_atoms": heavy,
                        "per_atom": round(best / heavy, 4),
                        "all_poses": [round(float(x), 3) for x in e[:, 0]]}

    (HERE / "scores.json").write_text(json.dumps(rows, indent=1))
    print(f"{'ligand':24s} {'score':>8s} {'heavy':>6s} {'per atom':>9s}")
    for n, r in sorted(rows.items(), key=lambda t: t[1]["score"]):
        print(f"{n:24s} {r['score']:8.2f} {r['heavy_atoms']:6d} {r['per_atom']:9.3f}")
    print("\nrank by raw score:      " + ", ".join(
        n for n, _ in sorted(rows.items(), key=lambda t: t[1]["score"])))
    print("rank by score per atom: " + ", ".join(
        n for n, _ in sorted(rows.items(), key=lambda t: t[1]["per_atom"])))


if __name__ == "__main__":
    main()
