#!/bin/bash
# Two conda environments. Docking is CPU, the proton affinity wants a GPU.
set -euo pipefail

conda create -y -n dock python=3.11
conda install -y -n dock -c conda-forge openbabel numpy rdkit openmm pdbfixer
# gemmi is a meeko dependency that meeko 0.8.0 does not declare, so pip will not pull it
# in and the first import fails. Install it explicitly.
conda run -n dock pip install vina meeko gemmi

# GPU4PySCF: match the wheel to the CUDA version on your node (nvidia-smi).
conda create -y -n qmgpu python=3.11
conda run -n qmgpu pip install pyscf gpu4pyscf-cuda12x cupy-cuda12x

echo
echo "built. point run.sbatch at them if they are not in ~/miniforge3/envs:"
echo "  DOCK_ENV=... DFT_ENV=... sbatch run.sbatch"

# DeepCCS is a 2019 model: python 3.7 and TensorFlow 1.15, so it needs its own environment and
# cannot share either of the others. It is CPU only, which is fine, because a prediction takes
# seconds and there is no GPU build of TF 1.15 worth the trouble.
# Do not pin numpy here. Pinning it alongside TF 1.15 makes the solve impossible.
conda create -y -n deepccs python=3.7
conda install -y -n deepccs -c conda-forge tensorflow=1.15 keras=2.3.1 h5py'<3' pandas
conda run -n deepccs pip install 'protobuf<3.21'
