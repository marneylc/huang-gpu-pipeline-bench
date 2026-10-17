#!/bin/bash -l
# The same docking, on 1,997 compounds instead of seven, with Uni-Dock on however many GPUs
# the job is given.
#
# This is the part that actually needs a GPU. Seven ligands do not: run_local.sh does those on
# CPU in minutes. The GPU case is throughput, and throughput is what a screening library needs.
# Every ligand is independent, so the work splits across GPUs with no communication between
# them, and the wall time should fall close to linearly with GPU count.
#
#   sbatch --gres=gpu:1 scale_gpu.sbatch     then 2, then 4, and compare ligands_per_s
#SBATCH -J acetyltropine_scale
#SBATCH -p preempt-gpu.q
#SBATCH -A YOUR_ACCOUNT
#SBATCH --gres=gpu:2
# Tesla K40/K80 nodes: driver too old for the CUDA 12 builds. Cluster-specific.
#SBATCH -x gpu001,bosun01
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -o scale_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
DOCK_ENV="${DOCK_ENV:-$HOME/miniforge3/envs/dock}"

[ -f receptor.pdbqt ] || { echo "run prep.py first"; exit 1; }
"$DOCK_ENV"/bin/python prep_library.py

# Which GPUs did Slurm actually give us? nvidia-smi -L lists every card on the node, which is
# not the same thing: a gpu:2 job on a four-GPU node would see four and try to use cards it was
# not allocated. Slurm puts the allocated ones in CUDA_VISIBLE_DEVICES, so take them from there.
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
  IFS=',' read -ra GPUS <<< "$CUDA_VISIBLE_DEVICES"
else
  mapfile -t GPUS < <(nvidia-smi -L | awk '{print NR-1}')
fi
NGPU=${#GPUS[@]}
echo "node $(hostname)  arch $(uname -m)  GPUs allocated: $NGPU (${GPUS[*]})"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader

read -r CX CY CZ < box.txt
ls library/*.pdbqt > all.txt
split -n l/"$NGPU" -d all.txt part_
t0=$(date +%s)
for ((g=0; g<NGPU; g++)); do
  mkdir -p out_g$g
  CUDA_VISIBLE_DEVICES=${GPUS[g]} "$DOCK_ENV"/bin/unidock \
    --receptor receptor.pdbqt --ligand_index part_0$g --dir out_g$g \
    --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
    --size_x 22 --size_y 22 --size_z 22 \
    --search_mode balance --scoring vina --num_modes 1 --seed 1 > dock_g$g.log 2>&1 &
done
wait
n=$(cat out_g*/*.pdbqt 2>/dev/null | grep -c "VINA RESULT" || true)
dt=$(( $(date +%s) - t0 ))
echo "DOCKED  ligands=$n  gpus=$NGPU  wall_s=$dt  ligands_per_s=$(awk -v n=$n -v t=$dt 'BEGIN{printf "%.1f", n/t}')"
echo
echo "Now where does acetyltropine sit among them?"
"$DOCK_ENV"/bin/python rank.py
