#!/bin/bash
# Two conda environments. Docking is CPU, the proton affinity wants a GPU.
set -euo pipefail

conda create -y -n dock python=3.11
conda install -y -n dock -c conda-forge openbabel numpy rdkit openmm pdbfixer
conda run -n dock pip install vina meeko

# GPU4PySCF: match the wheel to the CUDA version on your node (nvidia-smi).
conda create -y -n qmgpu python=3.11
conda run -n qmgpu pip install pyscf gpu4pyscf-cuda12x cupy-cuda12x

echo
echo "built. point run.sbatch at them if they are not in ~/miniforge3/envs:"
echo "  DOCK_ENV=... DFT_ENV=... sbatch run.sbatch"
