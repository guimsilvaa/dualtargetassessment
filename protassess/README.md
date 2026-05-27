# ProtAssess

ProtAssess is a python tool designed for **ensemble cavity analysis and comparison of protein binding pockets** using experimental protein structures. It makes the very good use of `pyKVFinder` cavity detection engine and it was developed to compare and analyze:

* ensembles of holo protein structures from distinct target families
* apo vs holo conformational changes
* pocket similarities/differences
* cavity physicochemical profiles

ProtAssess automatically performs:

1. Protein cleaning and preprocessing
2. Cavity detection
3. Ligand-proximal cavity selection
4. Cavity descriptor extraction
5. Statistical comparison
6. Multivariate analysis
7. Cavity visualization generation

---

## 1. Protein cleaning

Each structure is automatically cleaned using BioPython:

* removal of waters
* removal of ions/metals
* removal of co-crystallized ligands
* removal of alternate conformers
* removal of hydrogens
* extraction of a single biologically relevant monomer chain

### Ligand-proximal chain/monomer selection

~ reads a reference ligand SDF (within each target folder)
~ computes the ligand centroid
~ identifies the protein chain closest to the ligand
~ keeps only that chain

---

## 2. Cavity detection

Cavities are detected using `pyKVFinder` and computing:

* cavity volume
* depth
* hydropathy
* surface area
* cavity residues
* voxelized cavity geometries

---

## 3. Ligand-proximal cavity selection

`pyKVFinder` detects multiple cavities per structure. ProtAssess automatically:

* reconstructs cavity voxel coordinates
* computes cavity centroids
* measures distances to the ligand centroid
* selects the cavity closest to the ligand

---

## 4. Descriptor extraction

ProtAssess extracts geometric, topological and physicochemical cavity descriptors.

| Descriptor           | Meaning                          |
| -------------------- | -------------------------------- |
| volume               | cavity volume (Å³)               |
| surface_area         | cavity surface area              |
| avg_depth            | average cavity depth             |
| max_depth            | maximum cavity depth             |
| avg_hydropathy       | average cavity hydropathy        |
| n_interface_residues | number of cavity-lining residues |
| surface_volume_ratio | cavity exposure metric           |
| compactness          | buriedness/compactness metric    |
| sphericity           | globularity metric               |
| hydrophobic_fraction | fraction of hydrophobic residues |
| aromatic_fraction    | fraction of aromatic residues    |
| charged_fraction     | fraction of charged residues     |

Detailed descriptor definitions and formulas are available in: `descriptors_info.md`

---

## 5. Statistical analysis

* mean and standard deviation
* descriptor distribution through violin plots
* performs Welch t-tests p-values (Welch’s t-test is used instead of the standard Student t-test because it does not assume equal variances or equal sample sizes between ensembles, making it more appropriate for comparative structural datasets derived from heterogeneous protein ensembles)
* effect size ranking using Cohen’s d statistic (Cohen’s d effect size quantifies the magnitude of separation between target ensembles, enabling ranking descriptors according to their discriminatory power)

---

## 6. Multivariate Analysis

* hierarchical clustering dendrograms (separate structures according to overall cavity descriptor similarity)
* PCA (projects multidimensional descriptor space into two dimensions: scatter plot, feature loadings and explained variance)
* descriptor heatmaps (generated using z-score normalized descriptors and hierarchical row clustering)
* radar plots (summarize the mean z-score normalized descriptor profiles for each target ensemble)
* descriptor importance (ranks descriptors according to absolute Cohen’s d effect size)

---

## 7. Cavity visualization

ProtAssess exports cavity voxels as pseudoatom PDB files.

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

# Citation

If you use ProtAssess, please cite:

* pyKVFinder

and reference this ProtAssess repository.
