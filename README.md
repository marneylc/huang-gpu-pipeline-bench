# GPU pipeline benchmark for the Huang Complex seed fund application

This benchmarks the three GPU stages of a docking-plus-mass-spectrometry prediction pipeline: docking (Uni-Dock), a machine-learned potential (AIMNet2), and DFT (GPU4PySCF). It runs on Novus x86 GPU nodes (A40, RTX) and on the ARM GH200 Grace-Hopper node. Every molecule is independent, so each GPU gets its own share of the work, and throughput should scale close to linearly with GPU count.

Luke Marney, Oregon State University, October 2026. Test code for the 2026 OSU Huang Complex Supercomputing Seed Fund application "From target to detector, predicting the ions that therapeutic small molecules actually form in mass spectrometry". It was developed on the OSU Novus cluster.

## How the pipeline uses multiple GPUs

Every docking pose, conformer and DFT calculation is independent, so the pipeline uses many GPUs by running one worker per GPU on its own share of the work, with no communication between GPUs. `novus_gpu_benchmark.sbatch` shows the pattern for all three stages.

1. It counts the GPUs Slurm gave the job (`nvidia-smi -L`).
2. It splits the input list into that many parts (`split -n l/$NGPU`).
3. It launches one background process per GPU, each pinned with `CUDA_VISIBLE_DEVICES=$g`, and waits for all of them.

On the Huang Complex the same pattern scales from one node to many by giving each Slurm array task or node its own slice of the list. A job layer that packs calculations onto free GPUs, checkpoints results and retries failed convergence is part of the proposed work.

## Files

| File | What it does |
|---|---|
| `setup_envs.sh` | Builds the conda environments. On x86 they come from prebuilt packages. On ARM, GPU4PySCF and Uni-Dock are built from source. |
| `make_inputs.sbatch` | CPU job. Downloads the 5-HT2A receptor (PDB 6A93), prepares the receptor and 1,997 ligands, and embeds ML conformers. |
| `novus_gpu_benchmark.sbatch` | The benchmark. Runs all three stages on however many GPUs the job gets and reports the timings. |
| `dock/prep.py` | Receptor and ligand preparation (PDBFixer, Open Babel, RDKit, Meeko). |
| `dock/redock.py`, `dock/rmsd.py` | Redocking check for risperidone, and RMSD to the crystal pose. |
| `ml/ml_bench.py` | AIMNet2 conformer embedding and energy-plus-force throughput. |
| `sp.py`, `withaferin_a_mh.xyz`, `withaferin_a.xyz` | GPU4PySCF single points on withaferin A, the main withanolide of ashwagandha (*Withania somnifera*), as the neutral molecule (72 atoms) and protonated at the C-1 ketone (73 atoms). The two energies give a proton affinity, the kind of quantity the pipeline uses to rank ion forms. |
| `tested_ligands_5HT2A.csv` | Uni-Dock scores for the docked ligands. |
| `ligands_2000.smi` | 2,000 drug-like compounds (MW 180 to 500) from public MS/MS libraries. |

## Running on Novus

These notes are specific to the OSU Novus cluster. Replace `YOUR_ACCOUNT` with your Slurm account.

- Jobs must be submitted from `/scratch/$USER`. Submitting from home is rejected. Scratch is purged after 90 days, so the setup script keeps source builds in `~/src`.
- The login node has no GPU driver, so GPU code can only be tested inside a job.
- Exclude `gpu001` (Tesla K40, driver 470). Its driver is too old for CUDA 12.
- Point `TMPDIR` and `PIP_CACHE_DIR` at scratch, because `/tmp` is a small tmpfs. The scripts already do this.

## Steps

```bash
mkdir -p /scratch/$USER/pipeline_bench && cp -r novus_pipeline_bench/* /scratch/$USER/pipeline_bench/
cd /scratch/$USER/pipeline_bench

# 1. Build the software (about 20 min on x86, about 1 hour on ARM)
srun -p preempt-gpu.q -A YOUR_ACCOUNT -x gpu001 --gres=gpu:1 -c 16 --mem=32G -t 2:00:00 bash setup_envs.sh   # x86
srun -p cos-arm.q     -A YOUR_ACCOUNT           --gres=gpu:1 -c 16 --mem=64G -t 3:00:00 bash setup_envs.sh   # ARM GH200

# 2. Make the inputs (needs the x86 environments)
sbatch make_inputs.sbatch

# 3. Run the benchmark
sbatch novus_gpu_benchmark.sbatch                        # 2x A40 by default
sbatch --gres=gpu:4 novus_gpu_benchmark.sbatch           # more GPUs, if a node has them
sbatch -p cos-arm.q --gres=gpu:1 novus_gpu_benchmark.sbatch   # ARM GH200
```

## ARM build notes (aarch64)

None of the three packages installs cleanly from prebuilt packages on aarch64.

- **GPU4PySCF**: there is no aarch64 wheel. Build v1.8.1 from source with conda-forge CUDA 12.6 and `CUDA_ARCHITECTURES=90`, then `pip install pyscf-dispersion`.
- **Uni-Dock**: there is no aarch64 conda package. Build from GitHub with conda-forge CUDA, Boost, and `libcurand-dev`. Without `libcurand-dev`, the build fails on `curand_kernel.h`.
- **AIMNet2**: the PyTorch cu126 aarch64 wheel works. Its dependency NVIDIA Warp needs two fixes. Replace the pip `warp-lang` (it needs glibc 2.35, and the node has 2.34) with the conda-forge build. Then set `warp.config.use_precompiled_headers = False`, because torch's NVRTC 12.6 rejects `-pch`. `ml_bench.py` already sets this.

## Reference results (October 2026)

| Stage | 1x A40 | 2x A40 | 1x GH200 |
|---|---|---|---|
| Docking, 1,997 ligands | 145 s | 84 s | 89 s |
| AIMNet2, 9,980 conformers | 10.2 s | 5.4 s | 13.2 s |
| DFT, 51-atom cation / radical (see note) | 85 s / 249 s | 2.01x on 4 calcs | 3.7 s / 10.2 s |

The DFT timings in this row were measured on a 51-atom organic cation and its radical that are not distributed here. The withaferin A test molecule is larger (72 and 73 atoms), so its calculations take longer. For comparison, ORCA on 16 CPU cores takes 257 s / 428 s for the same DFT. GH200 and A40 energies are identical. Docking scores agree between architectures to 0.12 kcal/mol on average.
