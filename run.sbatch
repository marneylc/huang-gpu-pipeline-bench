#!/bin/bash -l
# The whole example, start to finish. Submit from a scratch directory:
#   sbatch run.sbatch
#
# Docking is CPU here (AutoDock Vina). The proton affinity is the GPU part, and the two single
# points are independent, so they go on separate GPUs at the same time. That is the same
# pattern scale_gpu.sbatch uses on thousands of ligands, just with two pieces of work instead.
# Ask for one GPU and it still runs; they simply take turns.
#SBATCH -J acetyltropine
#SBATCH -p preempt-gpu.q
#SBATCH -A YOUR_ACCOUNT
#SBATCH --gres=gpu:2
# gpu001 and bosun01 are Tesla K40/K80: their driver is too old for the CUDA 12 build
# of GPU4PySCF, and a job that lands there dies on import. Cluster-specific; delete
# this line elsewhere.
#SBATCH -x gpu001,bosun01
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
# Which GPUs did Slurm actually give us? nvidia-smi -L lists every card on the node, which is
# not the same thing: a gpu:2 job on a four-GPU node would see four and try to use cards it was
# not allocated. Slurm puts the allocated ones in CUDA_VISIBLE_DEVICES, so take them from there.
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
  IFS=',' read -ra GPUS <<< "$CUDA_VISIBLE_DEVICES"
else
  mapfile -t GPUS < <(nvidia-smi -L | awk '{print NR-1}')
fi
NGPU=${#GPUS[@]}
echo "GPUs allocated: $NGPU (${GPUS[*]})"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader
t0=$(date +%s)
# the cation and the neutral do not depend on each other, so run them at once when there is
# more than one GPU, and one after the other when there is not
CUDA_VISIBLE_DEVICES=${GPUS[0]} "$DFT_ENV"/bin/python sp.py acetyltropine_mh.xyz 1 1 cation \
  > dft_cation.txt 2>&1 &
CUDA_VISIBLE_DEVICES=${GPUS[1 % NGPU]} "$DFT_ENV"/bin/python sp.py acetyltropine.xyz 0 1 neutral \
  > dft_neutral.txt 2>&1 &
wait
cat dft_cation.txt dft_neutral.txt > dft_energies.txt
cat dft_energies.txt
echo "DFT  calcs=2  gpus=$NGPU  wall_s=$(( $(date +%s) - t0 ))"
"$DFT_ENV"/bin/python proton_affinity.py dft_energies.txt

echo; echo "=== 5. how it breaks in the mass spectrometer"
"$DOCK_ENV"/bin/python fragment.py

echo; echo "=== 6. what else arrives at the same mass"
"$DOCK_ENV"/bin/python interference.py
