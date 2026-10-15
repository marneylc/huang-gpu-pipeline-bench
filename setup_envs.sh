#!/bin/bash -l
# Build the three GPU environments used by novus_gpu_benchmark.sbatch.
# Tested on OSU Novus, 2026-10-03: x86_64 A40 nodes (driver 560, CUDA 12.6)
# and the aarch64 GH200 node cos-gh01 (glibc 2.34).
#
# Builds must run on a node of the target architecture, so the ARM build has
# to run on the GH200 itself. From /scratch/$USER:
#   x86:  srun -p preempt-gpu.q -A YOUR_ACCOUNT -x gpu001 --gres=gpu:1 -c 16 --mem=32G -t 2:00:00 bash setup_envs.sh
#   ARM:  srun -p cos-arm.q     -A YOUR_ACCOUNT           --gres=gpu:1 -c 16 --mem=64G -t 3:00:00 bash setup_envs.sh
# Each step is skipped if its environment already exists, so re-running is safe.

set -euo pipefail
ARCH=$(uname -m)
export TMPDIR=/scratch/$USER/tmp PIP_CACHE_DIR=/scratch/$USER/pipcache_$ARCH   # login/compute /tmp is a small tmpfs
mkdir -p "$TMPDIR" "$PIP_CACHE_DIR"
SRC=$HOME/src; mkdir -p "$SRC"   # source builds live in home: /scratch is purged after 90 days

if [ "$ARCH" = x86_64 ]; then MF=${MF:-~/miniforge3}; else MF=${MF:-~/miniforge3-arm}; fi   # override with MF=/path to install elsewhere
if [ ! -x "$MF/bin/mamba" ]; then
  curl -sL -o "$TMPDIR/mf.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$ARCH.sh"
  bash "$TMPDIR/mf.sh" -b -p "$MF"
fi
have() { [ -d "$MF/envs/$1" ]; }
CUDA_LIBS="cupy-cuda12x[ctk] nvidia-nvjitlink-cu12 nvidia-cublas-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12 nvidia-curand-cu12 nvidia-cufft-cu12 nvidia-cuda-runtime-cu12 nvidia-cuda-nvrtc-cu12"

########## x86_64: prebuilt packages ##########
if [ "$ARCH" = x86_64 ]; then
  # DFT: GPU4PySCF wheels exist for x86. CuPy needs the CUDA runtime libs from pip,
  # because there is no system CUDA toolkit on Novus.
  have qmgpu || { "$MF/bin/mamba" create -y -q -n qmgpu python=3.11
                  "$MF/envs/qmgpu/bin/pip" install -q gpu4pyscf-cuda12x cutensor-cu12 pyscf rdkit $CUDA_LIBS; }
  # Docking: Uni-Dock plus ligand/receptor preparation tools, all from conda-forge.
  have dock  || "$MF/bin/mamba" create -y -q -n dock -c conda-forge python=3.11 \
                  unidock meeko rdkit openbabel pdbfixer gemmi scipy
  # ML potential: PyTorch (CUDA 12.6 wheels to match driver 560) + AIMNet2.
  have mlpot || { "$MF/bin/mamba" create -y -q -n mlpot python=3.11
                  "$MF/envs/mlpot/bin/pip" install -q torch --index-url https://download.pytorch.org/whl/cu126
                  "$MF/envs/mlpot/bin/pip" install -q aimnet rdkit ase; }
fi

########## aarch64 (GH200): GPU4PySCF and Uni-Dock built from source ##########
if [ "$ARCH" = aarch64 ]; then
  # DFT: there is no aarch64 gpu4pyscf wheel, so build v1.8.1 for Hopper (sm_90).
  if ! have g4p; then
    "$MF/bin/mamba" create -y -q -n g4p -c conda-forge python=3.11 pyscf "cupy>=13" cuda-version=12.6 \
        cuda-nvcc cuda-cudart-dev libcublas-dev libcusolver-dev libcusparse-dev cutensor cmake make gxx_linux-aarch64 git
    rm -rf "$SRC/gpu4pyscf"
    git clone -q --depth 1 --branch v1.8.1 https://github.com/pyscf/gpu4pyscf.git "$SRC/gpu4pyscf"
    ( source "$MF/bin/activate" g4p && cd "$SRC/gpu4pyscf" &&
      cmake -B build -S gpu4pyscf/lib -DCUDA_ARCHITECTURES=90 -DCMAKE_CUDA_ARCHITECTURES=90 -DBUILD_LIBXC=ON &&
      cmake --build build -j "${SLURM_CPUS_PER_TASK:-16}" &&
      pip install -q --no-build-isolation --no-deps -e . &&
      pip install -q pyscf-dispersion )
  fi
  # Docking: there is no aarch64 unidock on conda-forge, so build from GitHub.
  # libcurand-dev is required (curand_kernel.h).
  if ! have udock; then
    "$MF/bin/mamba" create -y -q -n udock -c conda-forge cuda-version=12.6 cuda-nvcc cuda-cudart-dev \
        libcurand-dev cmake make gxx_linux-aarch64 libboost-devel git
    rm -rf "$SRC/Uni-Dock"
    git clone -q --depth 1 https://github.com/dptech-corp/Uni-Dock.git "$SRC/Uni-Dock"
    ( source "$MF/bin/activate" udock && cd "$SRC/Uni-Dock/unidock" &&
      cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=90 &&
      cmake --build build -j "${SLURM_CPUS_PER_TASK:-16}" &&
      ln -sf "$SRC/Uni-Dock/unidock/build/unidock" "$MF/envs/udock/bin/unidock" )
  fi
  # ML potential: the aarch64 PyTorch cu126 wheel works as-is. AIMNet2's dependency
  # NVIDIA Warp needs two fixes:
  #   - the pip warp-lang wheel needs glibc 2.35 (node has 2.34), so use the conda-forge build
  #   - set warp.config.use_precompiled_headers = False at runtime, because torch's NVRTC 12.6
  #     rejects -pch (ml_bench.py does this)
  if ! have mlpot; then
    "$MF/bin/mamba" create -y -q -n mlpot python=3.11
    "$MF/envs/mlpot/bin/pip" install -q torch --index-url https://download.pytorch.org/whl/cu126
    "$MF/envs/mlpot/bin/pip" install -q aimnet rdkit ase
    "$MF/envs/mlpot/bin/pip" uninstall -y -q warp-lang
    "$MF/bin/mamba" install -y -q -n mlpot -c conda-forge warp-lang
  fi
fi

echo "Environments in $MF/envs:"; ls "$MF/envs"
