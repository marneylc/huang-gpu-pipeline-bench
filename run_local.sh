#!/bin/bash
# Everything except the proton affinity, on one machine, no scheduler.
#
# Docking is CPU only, so this does not need a GPU at all. The DFT step is the only part that
# wants one, and it is left out here; run it on a cluster with run.sbatch, or drop it.
#
#   bash run_local.sh [n_cpus]
set -euo pipefail
cd "$(dirname "$0")"
NCPU="${1:-0}"          # 0 lets Vina pick

python prep.py
echo
python dock.py "$NCPU"
echo
python redock.py
echo
python protonation.py
echo
python interference.py
echo
echo "fragmentation needs the CFM-ID container; it fetches it the first time:"
echo "  python fragment.py"
echo
echo "the proton affinity step needs a GPU:"
echo "  python sp.py acetyltropine_mh.xyz 1 1 cation"
echo "  python sp.py acetyltropine.xyz    0 1 neutral"
