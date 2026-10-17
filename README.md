# Acetyltropine, from a target to a detector

![From a target to a detector: acetyltropine from ashwagandha root, with every stage of the pipeline predicted before a sample is run](docs/pipeline_slide.png)

A small, complete worked example of the question this project asks: **before you spend a year trying to measure a molecule, can you work out whether it is measurable at all?**

The compound is acetyltropine, a tropane alkaloid of ashwagandha root (*Withania somnifera*). It is a good example because it is small, it is real, and every stage gives an answer in minutes rather than days.

The code here covers the first two stages on that slide, the two that need a cluster:

1. **Does it bind?** Dock it into the M1 muscarinic receptor, the switch that atropine and scopolamine block.
2. **Does it take a charge?** Compute its proton affinity, which is what decides whether a MALDI matrix can ionize it.

There is a third thing in between, and it is the part worth your attention:

3. **Does the method work at all?** Redock a ligand whose real position is already known, and see whether the software finds it.

## Run it

On a cluster, all of it including the proton affinity:

```
bash setup_env.sh          # two conda environments, once
sbatch run.sbatch          # edit the account line first
```

On one machine with no scheduler, everything except the proton affinity. Docking is CPU only, so this does not need a GPU:

```
bash setup_env.sh          # the `dock` environment is the only one you need
conda activate dock
bash run_local.sh          # about ten minutes on 16 cores
```

Each step also runs on its own:

```
python prep.py             # fetch 5CXV, build the receptor and seven ligands
python dock.py             # dock them, report two kinds of score
python redock.py           # check the setup against a pose we already know
python protonation.py      # what getting the protonation wrong costs
python sp.py acetyltropine_mh.xyz 1 1 cation    # needs a GPU
```

## The three things this example is really about

**A docking score is not evidence that the pose is right.** `redock.py` takes tiotropium, whose position in 5CXV was measured by crystallography, and docks it back in. It lands 0.97 angstrom from where it actually sits, which passes. If it had not, nothing else the docking says could be trusted, and the scores would have looked perfectly reasonable anyway.

**Getting the protonation wrong is invisible in the score.** Every tropane here except tiotropium is a tertiary amine with a pKa near 9.5, so at body pH it carries a proton, and that ammonium makes the salt bridge that anchors it. `prep.py` builds them that way. `protonation.py` builds them both ways, changes nothing else, and docks both:

| ligand | protonated | neutral | score moves | pose moves |
|---|---|---|---|---|
| cocaine | -9.40 | -8.93 | -0.47 | **3.22 A** |
| 3-tigloyloxytropane | -7.72 | -7.79 | +0.07 | **2.37 A** |
| atropine | -9.14 | -9.29 | +0.16 | 1.49 A |
| scopolamine | -9.45 | -9.33 | -0.12 | 1.37 A |
| tropine | -5.72 | -5.98 | +0.26 | 0.69 A |
| acetyltropine | -6.51 | -6.50 | -0.01 | 0.08 A |
| tiotropium | -9.29 | -9.29 | +0.00 | no basic N |

The score never moves half a kcal/mol. The pose moves anywhere from nothing to three angstrom, and the two columns have nothing to do with each other. Acetyltropine barely moves because it is small and rigid; the flexible esters swing right round. So the answer is not "protonation always matters", it is that the error is invisible in the score and you cannot guess which case you are in. You check.

**How you keep score decides who wins.** `dock.py` prints the raw score and the score per heavy atom, and prints both rankings. They come out exactly reversed:

```
rank by raw score:      scopolamine, tiotropium, atropine, 3-tigloyloxytropane, acetyltropine, tropine
rank by score per atom: tropine, acetyltropine, 3-tigloyloxytropane, atropine, scopolamine, tiotropium
```

Bigger molecules make more contacts, so they climb the raw ranking almost automatically. Acetyltropine is 13 heavy atoms, about half the size of the drugs it is up against. Neither column is wrong. Reporting only one of them is.

## Files

| File | What it does |
|---|---|
| `prep.py` | Downloads PDB 5CXV, prepares chain A with PDBFixer and Open Babel, and builds the ligands from `ligands.smi` with their ring nitrogens protonated. |
| `ligands.smi` | Seven tropanes: acetyltropine, tropine and 3-tigloyloxytropane from the root, plus tiotropium, atropine, scopolamine and cocaine as known tropane ligands to measure against. |
| `dock.py` | AutoDock Vina, box centred on the crystal ligand. Writes `scores.json` and the poses. |
| `redock.py` | The correctness check described above. Symmetry-aware RMSD, computed in place without realignment, because realigning would quietly forgive a pose in the wrong place. |
| `protonation.py` | Docks every ligand protonated and neutral and tabulates what changed. |
| `sp.py`, `acetyltropine{,_mh}.xyz` | GPU4PySCF single points on the neutral molecule (30 atoms) and the ring-nitrogen protomer (31 atoms). The difference is the proton affinity. |
| `run.sbatch` | All of it, in one Slurm job. |
| `run_local.sh` | Everything except the proton affinity, on one machine, no scheduler. |

## What it needs

A machine with GPUs and Slurm. AutoDock Vina, RDKit, Meeko, PDBFixer, Open Babel and Open Babel's `obabel` on PATH for the docking side; PySCF and GPU4PySCF for the proton affinity. `setup_env.sh` builds both.

Developed on Oregon State's Novus cluster for the 2026 OSU Huang Complex Supercomputing Seed Fund application, "From target to detector, predicting the ions that therapeutic small molecules actually form in mass spectrometry". Luke Marney, Oregon State University.

The acetyltropine chemistry it rests on is in Dutta et al., *J Ethnopharmacol* **357**:120905 (2026), which synthesized the compound, separated it from its isomer acetyl exotropine by retention time and 2D-NMR, and showed the root extract's effect in *Drosophila*.
