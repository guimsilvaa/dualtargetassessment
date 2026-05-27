# ProtAssess

ProtAssess is a python tool designed for **ensemble cavity analysis and comparison of protein binding pockets** using experimental protein structures. It makes the very good use of `pyKVFinder` cavity detection engine and it was developed to compare and analyze:

* ensembles of holo protein structures from distinct target families
* apo vs holo conformational changes
* pocket similarities/differences
* cavity physicochemical profiles

ProtAssess automatically performs the following steps:

1. Protein cleaning and preprocessing
2. Cavity detection
3. Ligand-proximal cavity selection
4. Cavity descriptor extraction
5. Statistical comparison
6. Multivariate analysis
7. Cavity visualization generation

See more details about each step on [STEPS DETAILS](README_steps.md)

---

# Installation

## Option 1 environment:

```bash
conda env create -f environment.yml
conda activate protassess
```

## Option 2 manual conda install:

```bash
conda create -n protassess python=3.10
conda activate protassess
conda install -c conda-forge numpy pandas matplotlib scipy scikit-learn biopython tqdm rdkit
pip install pyKVFinder
```

## Option 3 uv:

```bash
uv venv protassess
source protassess/bin/activate
uv pip install numpy pandas matplotlib scipy scikit-learn biopython tqdm rdkit pyKVFinder
```

---

# Running ProtAssess

Inside the project directory:
```bash
python run_protassess.py
```
ProtAssess automatically detects the pair of folders beginning with *target_* and processes them as independent ensembles.

! please note it only supports 2 target folders per run!

---


# Folder Organization

Expected input structure:

```text
protassess/
├── target_folder1/
│   ├── structures.pdb
│   └── reference_ligand.sdf
│
├── target_folder2/
│   ├── structures.pdb
│   └── reference_ligand.sdf
│
├── run_protassess.py
├── detect_cavities.py
├── extract_descriptors.py
├── analyze_ensembles.py
├── cluster_pockets.py
└── utils.py
```

Each target folder:

* must contain a ensemble of aligned PDB structures
* must contain exactly one reference ligand SDF

---

# Output Structure

Each run automatically generates a folder *YYMMDD_target1_target2_protassess-out/* containing:

* cleaned structures
* cavity visualization files
* `results_summary.csv`
* violin plots for each descriptor
* descriptor tables/importances
* clustering plots (PCA, dendrogram, heatmap, radarplot)
* complete log file

---

# Citation

If you use ProtAssess, please cite:

* pyKVFinder

and reference this ProtAssess repository.
