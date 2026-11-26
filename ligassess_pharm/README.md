# LigAssess-Pharm

**Pharmacophore-perspective assessment of a target pair.**

`ligassess_pharm` derives a *consensus pharmacophore map* for a protein target from an ensemble of its co-crystallized ligands. Instead of describing one ligand, it asks where, in three-dimensional space, a target repeatedly demands a donor, an acceptor or an aromatic contact — and how strongly conserved each of those demands is.

Running it on two targets and comparing the resulting hotspot maps shows whether their recognition requirements can plausibly be satisfied by a single scaffold.

Part of [dualtargetassessment](../README.md).

---

## How it works

1. **Load ligands** — reads every `*.sdf` in the given folder, keeping only molecules with 3D conformers. Ligands are expected to be in a common reference frame (i.e. extracted from structurally aligned complexes).
2. **Extract features** — adds explicit hydrogens with coordinates and assigns pharmacophoric features using RDKit's `BaseFeatures.fdef`, mapped onto four canonical types: `Donor`, `Acceptor`, `Aromatic`, `Hydrophobe`.
3. **Cluster per feature type** — DBSCAN groups feature points of the same type that recur across ligands. Each feature type uses its own radius, since aromatic contacts tolerate more positional spread than directional hydrogen bonds:

   | Feature | `eps` (Å) | `min_samples` |
   |---|---|---|
   | Donor | 1.5 | 3 |
   | Acceptor | 1.5 | 3 |
   | Aromatic | 2.0 | 2 |
   | *(other)* | 2.0 | 2 |

4. **Score and rank** — each surviving cluster becomes a *hotspot*, ranked by the number of contributing feature points and assigned a normalized `size_score` (1.00 = most conserved hotspot in the map).
5. **Export** — hotspots are written as pseudoatoms so they can be opened in any molecular viewer alongside the ligands that generated them.

---

## Usage

```bash
python ligassess_pharm3.py PATH_TO_FOLDER
```

Where `PATH_TO_FOLDER` is a directory containing the aligned ligand SDF files for one target:

```bash
python ligassess_pharm3.py ache
python ligassess_pharm3.py gsk3b
```

> **Note:** output files are written to the **current working directory**, not to the input folder. Run the command from inside the target's directory, or move the outputs afterwards.

---

## Inputs

| Path | Contents |
|---|---|
| `all21ache.sdf` | 21 AChE co-crystallized ligands (aligned) |
| `all64gsk3b.sdf` | 64 GSK-3β co-crystallized ligands (aligned) |
| `ache/all21ache_uniprep.sdf` | prepared AChE ligand ensemble used for the published map |
| `gsk3b/all64gsk3b_uniprep.sdf` | prepared GSK-3β ligand ensemble used for the published map |
| `E20_A_604.sdf`, `EBW_A_601.sdf`, `H0L_A_601.sdf`, `H0R_A_601.sdf`, `H1R_A_601.sdf` | individual reference ligands |
| `BaseFeatures.fdef` | reference copy of the RDKit feature definitions |

The ligand ensembles were extracted from the SIENA binding site ensembles in [`../data/siena_ensembles.zip`](../data/siena_ensembles.zip).

---

## Outputs

| File | Description |
|---|---|
| `hotspots.csv` | All hotspots with feature type, centroid coordinates, number of contributing points, rank and `size_score` |
| `hotspots.sdf` | Hotspot centroids as pseudoatoms, for visualization |
| `ligands_plus_hotspots.sdf` | Input ligands plus the hotspot pseudoatoms in one file — open this to see the map in context |

Precomputed results for both targets are committed under `ache/` and `gsk3b/`.

By default only the top **2 hotspots per feature type** (`MAX_PER_FEATURE`) are exported to SDF, while `hotspots.csv` retains the complete ranked list.

---

## Tuning

Clustering behaviour is controlled by constants at the top of `ligassess_pharm3.py`:

```python
MIN_LIGANDS_PER_CLUSTER = 2   # minimum distinct ligands contributing to a hotspot
MAX_PER_FEATURE = 2           # hotspots per feature type exported to SDF
CLUSTER_PARAMS = {...}        # per-feature DBSCAN eps and min_samples
```

Increase `eps` for sparse or loosely aligned ensembles; raise `min_samples` to keep only strongly conserved hotspots.

---

## Dependencies

`numpy` · `pandas` · `scikit-learn` · `rdkit`

Covered by the environment in [`../protassess/environment.yml`](../protassess/environment.yml).
