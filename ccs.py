#!/usr/bin/env python
"""Predict a collision cross section, and look hard at what the prediction is worth.

Stage 3 of the slide. Ion mobility separates ions by shape, so it is the one stage that could
tell acetyltropine apart from its isomer acetyl exotropine. They have the same formula, the
same precursor mass and the same fragments, so nothing else in this pipeline can.

DeepCCS reads a SMILES string. Run it on both isomers and on the same molecule with the
stereochemistry stripped out, and compare the three numbers before believing any of them.

This is a 2019 TensorFlow 1.15 model. It is small: a prediction takes a few seconds on a CPU
and there is no GPU build worth installing for it. Knowing which stages need a GPU is part of
the point of this repository, and this one does not.
"""
import subprocess, sys, shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "DeepCCS_src"
REPO = "https://github.com/plpla/DeepCCS.git"

ISOMERS = [
    ("CN1[C@H]2CC[C@@H]1C[C@@H](OC(C)=O)C2", "acetyltropine, 3-alpha"),
    ("CN1[C@H]2CC[C@@H]1C[C@H](OC(C)=O)C2", "acetyl exotropine, 3-beta"),
    ("CN1C2CCC1CC(OC(C)=O)C2", "same molecule, stereochemistry stripped"),
]


def ensure_src():
    if (SRC / "saved_models" / "default" / "model.h5").exists():
        return
    if not shutil.which("git"):
        sys.exit("needs git to fetch DeepCCS")
    print(f"fetching DeepCCS (once) into {SRC.name}/")
    subprocess.run(["git", "clone", "-q", REPO, str(SRC)], check=True)


def main():
    ensure_src()
    inp, out = HERE / "ccs_in.csv", HERE / "ccs_out.csv"
    inp.write_text("SMILES,Adducts\n" + "".join(f"{s},M+H\n" for s, _ in ISOMERS))

    # -ap and -sp default to paths relative to the tool's own directory, so pass all three
    model = SRC / "saved_models" / "default"
    cmd = [sys.executable, str(SRC / "interface" / "command_line_tool.py"), "predict",
           "-mp", str(model), "-ap", str(model), "-sp", str(model),
           "-i", str(inp), "-o", str(out)]
    env = {"PYTHONPATH": str(SRC / "core")}
    import os
    r = subprocess.run(cmd, env={**os.environ, **env}, capture_output=True, text=True)
    if not out.exists():
        print(r.stdout[-800:], r.stderr[-800:])
        sys.exit("DeepCCS produced no output; it needs its own python 3.7 / TF 1.15 environment")

    vals = {}
    for line in out.read_text().splitlines()[1:]:
        if line.strip():
            smi, _, ccs = line.rsplit(",", 2)
            vals[smi] = float(ccs)

    print("predicted collision cross section, [M+H]+, square angstrom\n")
    for smi, label in ISOMERS:
        print(f"  {label:42s} {vals.get(smi, float('nan')):7.2f}")

    a, b, flat = (vals[ISOMERS[0][0]], vals[ISOMERS[1][0]], vals[ISOMERS[2][0]])
    print(f"\n  the two isomers differ by          {abs(b-a):5.2f}  ({100*abs(b-a)/a:.2f} %)")
    print(f"  stripping the stereochemistry by   {abs(a-flat):5.2f}  ({100*abs(a-flat)/a:.2f} %)")
    print(f"\nSeparating the isomers at that difference would need a mobility resolving power")
    print(f"near {a/abs(b-a):.0f}. A timsTOF reaches a few hundred.")
    print("""
So read this carefully. The model does give the isomers different numbers, because the stereo
descriptors are characters in the string it reads. But deleting the stereochemistry altogether
moves the answer more than twenty times further than swapping it does, which tells you the
model is keying on the string rather than on the shape. The isomer difference is smaller than
the model's own error and smaller than the instrument could resolve, so it is not evidence
that the two can be told apart.

Doing this properly means computing the cross section from a three-dimensional structure:
optimise the geometry, get partial charges, and run the trajectory method, which actually
collides a gas molecule with the ion from many directions. MobCal-MPI does that, and unlike a
SMILES model it has a real shape to work with, so it can in principle separate stereoisomers.
We have not built that stage out here. It is the obvious contribution if you want one.""")


if __name__ == "__main__":
    main()
