#!/usr/bin/env python
"""Turn the two single point energies into a proton affinity.

PA is the negative enthalpy change for B + H+ -> BH+, so with electronic energies it is
E(neutral) - E(cation), plus 5/2 RT for the translational enthalpy a bare proton carries.

This is the electronic value. A zero point correction would lower it by roughly 6 to 8
kcal/mol, which is systematic across molecules of this kind, so the ordering between compounds
survives without it but the absolute number runs high. Computing that correction needs
frequencies, which is a bigger calculation than this example is for.
"""
import re, sys
from pathlib import Path

HA_TO_KCAL = 627.5094740631
RT_TERM = 1.4816          # 5/2 RT at 298.15 K, the translational enthalpy of a bare proton


def read(tag, text):
    m = re.search(rf"^{tag}\s+\S+\s+natm=(\d+)\s+nao=\d+\s+E=(-?\d+\.\d+)\s+conv=(\w+)",
                  text, re.M)
    if not m:
        return None
    return int(m.group(1)), float(m.group(2)), m.group(3) == "True"


def main():
    # the energies are in whatever sp.py printed, so read the job log if one was given
    text = Path(sys.argv[1]).read_text() if len(sys.argv) > 1 else sys.stdin.read()
    neutral, cation = read("neutral", text), read("cation", text)
    if not neutral or not cation:
        sys.exit("could not find both 'neutral' and 'cation' lines from sp.py")
    for tag, (natm, e, conv) in (("neutral", neutral), ("cation", cation)):
        if not conv:
            sys.exit(f"{tag} did not converge; the energy is meaningless")
    pa = (neutral[1] - cation[1]) * HA_TO_KCAL + RT_TERM
    print(f"\nproton affinity of acetyltropine = {pa:.1f} kcal/mol  (electronic, no zero point)")
    print("A tertiary amine this basic sits above every common MALDI matrix, so ionizing it is")
    print("not the difficulty. Seeing it against the background is, which is stage 6.")


if __name__ == "__main__":
    main()
