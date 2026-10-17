#!/usr/bin/env python
"""Where does acetyltropine sit among the whole library, on each of the two scores?

This is the question the small example cannot answer. Seven ligands tell you nothing about
whether a score is good; two thousand give it a distribution to sit in.
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAMED = [l.split()[1] for l in (HERE / "ligands.smi").read_text().splitlines()
         if l.strip() and not l.startswith("#")]


def main():
    rows = {}
    for f in sorted(HERE.glob("out_g*/*.pdbqt")):
        txt = f.read_text()
        m = re.search(r"VINA RESULT:\s+(-?\d+\.\d+)", txt)
        if not m:
            continue
        heavy = sum(1 for l in txt.splitlines()
                    if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("H", "HD"))
        if heavy:
            rows[f.stem] = (float(m.group(1)), heavy)
    if not rows:
        return print("no docked poses found; run the docking step first")

    by_raw = sorted(rows, key=lambda k: rows[k][0])
    by_le = sorted(rows, key=lambda k: rows[k][0] / rows[k][1])
    n = len(rows)
    print(f"{n} ligands docked\n")
    print(f"{'ligand':24s} {'score':>7s} {'heavy':>6s} {'raw rank':>12s} {'per-atom rank':>15s}")
    for name in NAMED:
        if name not in rows:
            continue
        s, h = rows[name]
        print(f"{name:24s} {s:7.2f} {h:6d} {by_raw.index(name)+1:7d}/{n} "
              f"{by_le.index(name)+1:10d}/{n}")
    print("\nA small molecule can sit near the bottom on one and near the top on the other.")


if __name__ == "__main__":
    main()
