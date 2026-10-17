#!/usr/bin/env python
"""Predict how acetyltropine breaks in the mass spectrometer, with CFM-ID.

Stage 4 of the slide. Knowing the precursor mass is not enough for a targeted assay; you need
the transitions, which means knowing which fragments it actually gives.

CFM-ID ships as a container. This fetches it if Apptainer or Docker is available and tells you
what to do if neither is, rather than failing with an import error.
"""
import re, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SIF = HERE / "cfmid.sif"
IMAGE = "docker://wishartlab/cfmid:latest"
MODEL = "/trained_models_cfmid4.0/[M+H]+"
# CFM-ID reports three collision settings; the middle one is what a library search would use
SETTING = 1


def ensure_image():
    if SIF.exists():
        return "apptainer"
    if shutil.which("apptainer") or shutil.which("singularity"):
        exe = shutil.which("apptainer") or shutil.which("singularity")
        print(f"fetching {IMAGE} (once, about 2 GB)")
        subprocess.run([exe, "pull", str(SIF), IMAGE], check=True)
        return "apptainer"
    if shutil.which("docker"):
        return "docker"
    sys.exit("needs Apptainer or Docker:\n"
             f"  apptainer pull {SIF.name} {IMAGE}\n"
             "then run this again.")


def run(engine, smi):
    out = HERE / "cfmid_out"; out.mkdir(exist_ok=True)
    target = "/out/acetyltropine.txt"
    cmd = ([shutil.which("apptainer") or shutil.which("singularity"), "exec",
            "-B", f"{out}:/out", str(SIF)] if engine == "apptainer" else
           ["docker", "run", "--rm", "-v", f"{out}:/out", "wishartlab/cfmid:latest"])
    cmd += ["cfm-predict", smi, "0.001", f"{MODEL}/param_output.log",
            f"{MODEL}/param_config.txt", "1", target]
    subprocess.run(cmd, check=True)
    return out / "acetyltropine.txt"


def main():
    smi = next(l.split()[0] for l in (HERE / "ligands.smi").read_text().splitlines()
               if l.strip() and not l.startswith("#") and l.split()[1] == "acetyltropine")
    path = run(ensure_image(), smi)

    blocks, cur = {}, None
    for line in path.read_text().splitlines():
        if line.startswith("energy"):
            cur = int(line[6:]); blocks[cur] = []
        elif cur is not None and re.match(r"^\d", line):
            p = line.split()
            blocks[cur].append((float(p[0]), float(p[1])))
    peaks = sorted(blocks.get(SETTING, []), key=lambda t: -t[1])
    if not peaks:
        sys.exit("no peaks parsed from the CFM-ID output")

    known = {124.11: "loses acetic acid, the tropane core",
             142.12: "loses ketene, gives protonated tropine",
             184.13: "the whole molecule"}
    print(f"acetyltropine [M+H]+, CFM-ID middle collision setting\n")
    print(f"{'m/z':>10s} {'rel':>7s}   what it is")
    for mz, inten in peaks[:8]:
        note = next((v for k, v in known.items() if abs(mz - k) < 0.05), "")
        print(f"{mz:10.4f} {inten:7.1f}   {note}")
    print("\nThe two transitions a targeted assay would watch are 184.13 to 124.11 and\n"
          "184.13 to 142.12, and both come out of the prediction with no measurement.")


if __name__ == "__main__":
    main()
