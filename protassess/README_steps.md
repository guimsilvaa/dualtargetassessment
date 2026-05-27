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

Detailed descriptor definitions and formulas are available in [DESCRIPTORS DETAILS](README_descriptors.md)

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
