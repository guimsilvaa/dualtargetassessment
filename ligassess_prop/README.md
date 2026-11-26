# LigAssess-Prop

**Property-perspective assessment of a target pair.**

`ligassess_prop` compares the physicochemical and drug-likeness profiles of the ligand sets of two protein targets. Where the pocket and pharmacophore modules ask *what each target demands structurally*, this one asks a blunter question: **do the two targets' actives even occupy the same region of property space?**

When the distributions diverge consistently — in size, polarity, flexibility or sp³ character — a single drug-like scaffold satisfying both targets becomes correspondingly harder to find.

Part of [dualtargetassessment](../README.md).

---

## How it works

1. **Standardize** — every molecule is reduced to its parent fragment, normalized, uncharged and sanitized with RDKit's `MolStandardize`. Molecules that fail are dropped.
2. **Compute descriptors** — 15 properties per molecule (below).
3. **Compare distributions** — for each property and each target pair, a two-sided **Mann–Whitney U** test and a **Kolmogorov–Smirnov** test, plus **Cliff's delta** as a non-parametric effect size, classified into negligible / small / medium / large.
4. **Visualize** — a z-score normalized heatmap across all properties, plus kernel density plots generated **only for the properties that separate the two targets significantly**.

Non-parametric statistics are used throughout because ligand property distributions are typically skewed and the two sets differ greatly in size.

### Properties computed

`MW` · `logP` · `TPSA` · `HBD` · `HBA` · `RotBonds` · `FractionCSP3` · `HeavyAtoms` · `RingCount` · `AromaticRings` · `QED` · `FormalCharge` · `HeteroAtoms` · `LargestRingSize` · `FlexibilityIndex`

---

## Usage

```bash
python ligassess_prop.py pdb_prot1.sdf chembl_prot1.sdf pdb_prot2.sdf chembl_prot2.sdf
```

Exactly four SDF files are required, in this order: the **structural** (PDB-derived) and **bioactivity** (ChEMBL-derived) ligand sets of target 1, then the same two for target 2.

Reproducing the published comparison:

```bash
python ligassess_prop.py all21ache.sdf ache_2527actives.sdf all64gsk3b.sdf gsk3b_1646actives.sdf
```

> **Important:** target names are hardcoded at the top of `ligassess_prop.py` and must be edited when analyzing a different pair:
> ```python
> PROT1_NAME = "AChE"
> PROT2_NAME = "GSK3b"
> ```

---

## Inputs

| File | Target | Source | Contents |
|---|---|---|---|
| `all21ache.sdf` | AChE | PDB | 21 co-crystallized ligands |
| `ache_2527actives.sdf` | AChE | ChEMBL | 2,527 actives (IC₅₀ < 10 µM) |
| `all64gsk3b.sdf` | GSK-3β | PDB | 64 co-crystallized ligands |
| `gsk3b_1646actives.sdf` | GSK-3β | ChEMBL | 1,646 actives (IC₅₀ < 10 µM) |

---

## Outputs

| File | Description |
|---|---|
| `ligand_properties.csv` | Per-molecule property table, with target and dataset labels |
| `ligand_property_stats.csv` | Full statistics per property: medians, U and KS p-values, Cliff's delta, effect classification |
| `ligand_property_significant.csv` | Filtered to properties that separate the two targets significantly |
| `plots/<property>_density.png` | Kernel density distribution, both targets overlaid — written only for significantly separating properties |
| `plots/property_heatmap.png` | z-score normalized heatmap across all properties |

Precomputed outputs from the published run are committed in this directory and in `plots/`.

---

## Interpreting the results

Cliff's delta carries the signal worth reading first. Its **sign** says which target skews higher; its **magnitude** says whether the separation matters:

| \|delta\| | Interpretation |
|---|---|
| < 0.147 | negligible |
| 0.147 – 0.33 | small |
| 0.33 – 0.474 | medium |
| ≥ 0.474 | large |

Properties showing **medium or large** effects in **opposite directions** are the ones that constrain multitarget design: each target is pulling the scaffold toward a different region of property space.

---

## Dependencies

`numpy` · `pandas` · `matplotlib` · `seaborn` · `scipy` · `rdkit`

Covered by the environment in [`../protassess/environment.yml`](../protassess/environment.yml).
