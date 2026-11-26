# dualtargetassessment

Tools for assessing whether two protein targets can realistically be engaged by a **single drug-like scaffold**, by comparing them from three complementary perspectives:

| Perspective | Module | What it compares |
|---|---|---|
| **Pocket** | [`protassess/`](protassess/) | Binding site ensembles — cavity geometry, topology and physicochemistry |
| **Pharmacophore** | [`ligassess_pharm/`](ligassess_pharm/) | Spatial pharmacophoric hotspots derived from co-crystallized ligand ensembles |
| **Properties** | [`ligassess_prop/`](ligassess_prop/) | Physicochemical and drug-likeness property distributions of known actives |

---

## Associated publication

This repository contains the code and datasets supporting:

> Silva, G. M.; Klein, R.; Guerra, J. V. S.; Gastreich, M.; Alves, V. M.; Trivella, D. B. B.; Silva, C. H. T. P.; Muratov, E. N.; Tropsha, A.
> *Ultra-large virtual screening and target-pair assessment reveal conflicting requirements for single-scaffold drug-like GSK-3β/AChE dual inhibition.*
> Manuscript submitted. DOI to be added upon publication.

The case study throughout is the target pair **acetylcholinesterase (AChE)** and **glycogen synthase kinase-3β (GSK-3β)**, but all three modules are target-agnostic.

---

## Repository layout

```text
dualtargetassessment/
├── protassess/          # pocket-perspective assessment (see its own README)
├── ligassess_pharm/     # pharmacophore-perspective assessment
├── ligassess_prop/      # property-perspective assessment
├── data/
│   └── siena_ensembles.zip   # SIENA binding-site ensembles for AChE and GSK-3β
├── LICENSE
└── README.md
```

`data/siena_ensembles.zip` contains the two aligned binding-site ensembles (`siena_ache/`, `siena_gsk3b/`) used as input to `protassess`, each with its `ensemble/` PDB structures and the SIENA `alignment.txt`.

---

## Requirements

All modules run on Python 3.10+ and share a common scientific stack. The most complete environment file is `protassess/environment.yml`, which also covers the two `ligassess` modules:

```bash
conda env create -f protassess/environment.yml
conda activate protassess
```

Minimal per-module dependencies, if you prefer to install by hand:

| Module | Dependencies |
|---|---|
| `protassess` | numpy, pandas, matplotlib, scipy, scikit-learn, biopython, tqdm, rdkit, pyKVFinder |
| `ligassess_pharm` | numpy, pandas, scikit-learn, rdkit |
| `ligassess_prop` | numpy, pandas, matplotlib, seaborn, scipy, rdkit |

---

## Quick start

Each module is run independently from inside its own directory.

**1. Pocket perspective**

```bash
cd protassess
python run_protassess.py
```

Auto-detects exactly two folders named `target_*`, each holding an aligned PDB ensemble plus one reference ligand SDF. Full documentation in [`protassess/README.md`](protassess/README.md).

**2. Pharmacophore perspective**

```bash
cd ligassess_pharm
python ligassess_pharm3.py ache
```

Clusters RDKit pharmacophoric features across a ligand ensemble into consensus hotspots. See [`ligassess_pharm/README.md`](ligassess_pharm/README.md).

**3. Property perspective**

```bash
cd ligassess_prop
python ligassess_prop.py all21ache.sdf ache_2527actives.sdf all64gsk3b.sdf gsk3b_1646actives.sdf
```

Compares property distributions between the two targets' ligand sets. See [`ligassess_prop/README.md`](ligassess_prop/README.md).

---

## Reproducing the published analyses

The input files committed here are the exact ones used in the paper:

| Paper figure | Module | Inputs in this repository |
|---|---|---|
| Binding site ensemble comparison | `protassess` | `data/siena_ensembles.zip` |
| Pharmacophoric hotspot maps | `ligassess_pharm` | `all21ache.sdf`, `all64gsk3b.sdf` |
| Ligand property distributions | `ligassess_prop` | `all21ache.sdf`, `ache_2527actives.sdf`, `all64gsk3b.sdf`, `gsk3b_1646actives.sdf` |

Precomputed outputs are also included (`ligassess_pharm/*/hotspots.csv`, `ligassess_prop/ligand_propert*.csv`, `ligassess_prop/plots/`) so results can be inspected without rerunning anything.

---

## Citation

If you use these tools, please cite the associated publication above, together with the underlying methods:

* **pyKVFinder** — cavity detection engine used by `protassess`: https://github.com/LBC-LNBio/pyKVFinder
* **RDKit** — pharmacophoric feature definitions and descriptor calculation: https://www.rdkit.org
* **SIENA** — binding site ensemble generation: https://proteins.plus
