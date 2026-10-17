#!/bin/bash -l
# GPU benchmark for the target-to-detector pipeline on the OSU Novus cluster.
# Runs three stages on every GPU the job is given, one independent worker per GPU:
#   1. Docking      Uni-Dock, ligands split evenly across GPUs
#   2. ML potential AIMNet2 energies + forces on RDKit conformers
#   3. DFT          GPU4PySCF wB97X/def2-TZVP single points
# Throughput should scale close to linearly with GPU count, because every
# molecule is independent.
#
# Submit from scratch (Novus does not allow jobs from $HOME):
#   cd /scratch/$USER/huang_bench && sbatch novus_gpu_benchmark.sbatch
# Change GPU count with:  sbatch --gres=gpu:2 novus_gpu_benchmark.sbatch
#
# Build the software first with setup_envs.sh, and make the inputs with make_inputs.sbatch.
# Expected inputs in the working directory:
#   dock/receptor.pdbqt, dock/box.txt (center x y z), dock/all.txt (one ligand .pdbqt per line)
#   ml/ml_bench.py, ml/confs.npy
#   sp.py, acetyltropine_mh.xyz, acetyltropine.xyz

#SBATCH -J pipeline_gpu_bench
#SBATCH -p preempt-gpu.q          # A40 / RTX nodes; cos-arm.q for the GH200 (ARM)
#SBATCH -A YOUR_ACCOUNT    # your Slurm account
#SBATCH -x gpu001                 # Tesla K40 node: driver too old for CUDA 12
#SBATCH --gres=gpu:2
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -o pipeline_bench_%j.out

set -euo pipefail
export TMPDIR=/scratch/$USER/tmp; mkdir -p "$TMPDIR"

# Environments built by setup_envs.sh. The script picks the right set for the node's CPU.
if [ "$(uname -m)" = aarch64 ]; then
  DOCK_ENV=~/miniforge3-arm/envs/udock; ML_ENV=~/miniforge3-arm/envs/mlpot; DFT_ENV=~/miniforge3-arm/envs/g4p
else
  DOCK_ENV=~/miniforge3/envs/dock;      ML_ENV=~/miniforge3/envs/mlpot;      DFT_ENV=~/miniforge3/envs/qmgpu
fi

NGPU=$(nvidia-smi -L | wc -l)
echo "node $(hostname)  arch $(uname -m)  GPUs: $NGPU"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader

# ---------- 1. Docking ----------
cd dock
read -r CX CY CZ < box.txt
split -n l/"$NGPU" -d all.txt part_            # part_00, part_01, ...
t0=$(date +%s)
for ((g=0; g<NGPU; g++)); do
  mkdir -p out_g$g
  CUDA_VISIBLE_DEVICES=$g "$DOCK_ENV"/bin/unidock \
    --receptor receptor.pdbqt --ligand_index part_0$g --dir out_g$g \
    --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
    --size_x 22 --size_y 22 --size_z 22 \
    --search_mode balance --scoring vina --num_modes 1 --seed 1 > dock_g$g.log 2>&1 &
done
wait
nlig=$(cat out_g*/*.pdbqt 2>/dev/null | grep -c "VINA RESULT" || true)
dt=$(( $(date +%s) - t0 ))
echo "DOCKING  ligands=$nlig  wall_s=$dt  ligands_per_s=$(awk -v n=$nlig -v t=$dt 'BEGIN{printf "%.1f", n/t}')"
cd ..

# ---------- 2. ML potential ----------
cd ml
t0=$(date +%s)
for ((g=0; g<NGPU; g++)); do
  CUDA_VISIBLE_DEVICES=$g "$ML_ENV"/bin/python ml_bench.py run $g "$NGPU" 2>&1 | grep part &
done
wait
echo "ML  wall_s=$(( $(date +%s) - t0 ))  (includes model load)"
cd ..

# ---------- 3. DFT ----------
t0=$(date +%s)
for ((g=0; g<NGPU; g++)); do
  ( CUDA_VISIBLE_DEVICES=$g "$DFT_ENV"/bin/python sp.py acetyltropine_mh.xyz 1 1 gpu$g
    CUDA_VISIBLE_DEVICES=$g "$DFT_ENV"/bin/python sp.py acetyltropine.xyz 0 1 gpu$g ) &
done
wait
echo "DFT  calcs=$((2*NGPU))  wall_s=$(( $(date +%s) - t0 ))"
