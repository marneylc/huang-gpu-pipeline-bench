#!/bin/bash -l
# Build the benchmark inputs (CPU only). Run once from the benchmark directory on scratch,
# laid out as:  ./ligands_2000.smi  ./dock/prep.py  ./ml/ml_bench.py
#   sbatch make_inputs.sbatch
#SBATCH -J make_inputs
#SBATCH -p dri.q
#SBATCH -A YOUR_ACCOUNT    # your Slurm account
#SBATCH -c 32
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -o make_inputs_%j.out
set -euo pipefail
export TMPDIR=/scratch/$USER/tmp; mkdir -p "$TMPDIR"
MF=~/miniforge3   # x86 envs; the outputs are plain files usable on ARM too

# Receptor: serotonin 5-HT2A with risperidone (PDB 6A93). prep.py writes receptor.pdbqt,
# box.txt (risperidone centroid), and ligs/*.pdbqt from ../ligands_2000.smi via RDKit + Meeko.
cd dock
[ -f 6A93.pdb ] || curl -s -o 6A93.pdb https://files.rcsb.org/download/6A93.pdb
PATH=$MF/envs/dock/bin:$PATH python prep.py
[ -s receptor.pdbqt ] || $MF/envs/dock/bin/obabel recA_H.pdb -xr -O receptor.pdbqt
ls ligs/*.pdbqt > all.txt
echo "ligands: $(wc -l < all.txt)"
cd ..

# ML conformers: 10 RDKit conformers for each of the first 1,000 ligands -> ml/confs.npy
cd ml && $MF/envs/mlpot/bin/python ml_bench.py embed
