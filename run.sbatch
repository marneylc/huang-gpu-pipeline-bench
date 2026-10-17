#!/bin/bash -l
# The whole example, start to finish. Submit from a scratch directory:
#   sbatch run.sbatch
#
# Docking is CPU here (AutoDock Vina). Only the proton affinity needs a GPU, and on a molecule
# this size it is a couple of minutes, so the whole thing fits in one short job.
#SBATCH -J acetyltropine
#SBATCH -p preempt-gpu.q
#SBATCH -A YOUR_ACCOUNT
#SBATCH --gres=gpu:1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -o acetyltropine_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

DOCK_ENV="${DOCK_ENV:-$HOME/miniforge3/envs/dock}"
DFT_ENV="${DFT_ENV:-$HOME/miniforge3/envs/qmgpu}"

echo "=== 1. prepare the receptor and ligands"
"$DOCK_ENV"/bin/python prep.py

echo; echo "=== 2. dock them"
"$DOCK_ENV"/bin/python dock.py "$SLURM_CPUS_PER_TASK"

echo; echo "=== 3. does the setup reproduce a known pose?"
"$DOCK_ENV"/bin/python redock.py
echo; echo "=== 3b. what getting the protonation wrong costs"
"$DOCK_ENV"/bin/python protonation.py

echo; echo "=== 4. proton affinity of acetyltropine"
"$DFT_ENV"/bin/python sp.py acetyltropine_mh.xyz 1 1 cation
"$DFT_ENV"/bin/python sp.py acetyltropine.xyz    0 1 neutral
echo "proton affinity = (E_neutral - E_cation) * 627.5094 + 1.48 kcal/mol"
