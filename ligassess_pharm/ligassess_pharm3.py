#!/usr/bin/env python3

import os
import sys
import glob
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit.Geometry import Point3D

from sklearn.cluster import DBSCAN


# =========================
# PARAMETERS
# =========================
MIN_LIGANDS_PER_CLUSTER = 2
MAX_PER_FEATURE = 2


# Feature-specific clustering parameters
CLUSTER_PARAMS = {
    "Donor": {"eps": 1.5, "min_samples": 3},
    "Acceptor": {"eps": 1.5, "min_samples": 3},
    "Aromatic": {"eps": 2.0, "min_samples": 2},
}


# =========================
# FEATURE TYPE MAPPING
# =========================
FEATURE_MAP = {
    "Donor": ["Donor", "SingleAtomDonor"],
    "Acceptor": ["Acceptor", "SingleAtomAcceptor"],
    "Aromatic": ["Aromatic", "Arom4", "Arom5", "Arom6", "Arom7", "Arom8"],
    "Hydrophobe": ["Hydrophobe", "LumpedHydrophobe"]
}


def map_feature(ftype):
    for key, variants in FEATURE_MAP.items():
        if ftype in variants:
            return key
    return None


# =========================
# LOAD LIGANDS
# =========================
def load_ligands_from_folder(folder):
    sdf_files = glob.glob(os.path.join(folder, "*.sdf"))

    if not sdf_files:
        raise ValueError(f"No SDF files found in {folder}")

    mols = []

    for sdf in sdf_files:
        suppl = Chem.SDMolSupplier(sdf, removeHs=False)
        for mol in suppl:
            if mol is None:
                continue
            if not mol.GetNumConformers():
                continue
            mols.append(mol)

    return mols


# =========================
# FEATURE FACTORY
# =========================
def build_feature_factory():
    fdef = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')

    if not os.path.exists(fdef):
        raise FileNotFoundError("BaseFeatures.fdef not found.")

    return ChemicalFeatures.BuildFeatureFactory(fdef)


# =========================
# FEATURE EXTRACTION
# =========================
def extract_features(mol, factory, mol_id):
    mol = Chem.AddHs(mol, addCoords=True)

    feats = factory.GetFeaturesForMol(mol)
    data = []

    for f in feats:
        mapped_type = map_feature(f.GetType())
        if mapped_type is None:
            continue

        pos = f.GetPos()

        data.append({
            "type": mapped_type,
            "x": pos.x,
            "y": pos.y,
            "z": pos.z,
            "ligand_id": mol_id
        })

    return data


# =========================
# CLUSTERING
# =========================
def cluster_features(features, eps, min_samples):
    coords = np.array([[f["x"], f["y"], f["z"]] for f in features])

    if len(coords) == 0:
        return np.array([])

    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    return clustering.labels_


# =========================
# CLUSTER SUMMARY
# =========================
def summarize_clusters(features, labels):
    df = pd.DataFrame(features)
    df["cluster"] = labels

    clusters = []

    for c in set(labels):
        if c == -1:
            continue

        subset = df[df["cluster"] == c]

        ligand_ids = set(subset["ligand_id"])
        n_ligands = len(ligand_ids)

        if n_ligands < MIN_LIGANDS_PER_CLUSTER:
            continue

        centroid = subset[["x", "y", "z"]].mean().values

        clusters.append({
            "feature": subset["type"].iloc[0],
            "cluster_id": int(c),
            "count_points": int(len(subset)),
            "n_ligands": int(n_ligands),
            "x": float(centroid[0]),
            "y": float(centroid[1]),
            "z": float(centroid[2])
        })

    return pd.DataFrame(clusters)


# =========================
# CREATE HOTSPOT MOLECULE
# =========================
def create_hotspot_mol(cluster_df):
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(cluster_df))

    for i, row in cluster_df.iterrows():

        if row["feature"] == "Donor":
            atom = Chem.Atom(7)
        elif row["feature"] == "Acceptor":
            atom = Chem.Atom(8)
        else:
            atom = Chem.Atom(16)

        idx = mol.AddAtom(atom)
        conf.SetAtomPosition(idx, Point3D(row["x"], row["y"], row["z"]))

    mol.AddConformer(conf)
    mol = mol.GetMol()

    for i, row in cluster_df.iterrows():
        atom = mol.GetAtomWithIdx(i)

        atom.SetProp("Feature", row["feature"])
        atom.SetProp("CountPoints", str(row["count_points"]))
        atom.SetProp("NumLigands", str(row["n_ligands"]))
        atom.SetProp("ClusterID", str(row["cluster_id"]))
        atom.SetProp("Rank", str(row["rank"]))
        atom.SetProp("SizeScore", f"{row['size_score']:.2f}")

    mol.SetProp("_Name", "Pharmacophore_Hotspots")

    return mol


# =========================
# MAIN
# =========================
def main():
    if len(sys.argv) < 2:
        print("Usage: python ligassess_pharm.py PATH_TO_FOLDER")
        sys.exit(1)

    folder = " ".join(sys.argv[1:])

    print(f"[INFO] Loading ligands from: {folder}")
    mols = load_ligands_from_folder(folder)
    print(f"[INFO] Loaded {len(mols)} ligands")

    factory = build_feature_factory()

    print("[INFO] Extracting pharmacophore features...")
    all_features = []
    total_ligands = len(mols)

    for i, mol in enumerate(mols):
        feats = extract_features(mol, factory, i)
        all_features.extend(feats)

    print(f"[INFO] Total features collected: {len(all_features)}")

    # =========================
    # FEATURE-AWARE CLUSTERING
    # =========================
    all_clusters = []

    for ftype in set(f["type"] for f in all_features):
        subset = [f for f in all_features if f["type"] == ftype]

        params = CLUSTER_PARAMS.get(ftype, {"eps": 2.0, "min_samples": 2})

        labels = cluster_features(subset, params["eps"], params["min_samples"])

        cluster_df = summarize_clusters(subset, labels)

        if not cluster_df.empty:
            all_clusters.append(cluster_df)

    if not all_clusters:
        print("[WARNING] No clusters found.")
        sys.exit(0)

    cluster_df = pd.concat(all_clusters, ignore_index=True)

    # =========================
    # GLOBAL RANKING
    # =========================
    cluster_df = cluster_df.sort_values(
        by=["count_points"],
        ascending=False
    ).reset_index(drop=True)

    cluster_df["rank"] = cluster_df.index + 1
    cluster_df["size_score"] = cluster_df["count_points"] / cluster_df["count_points"].max()

    # Save full CSV
    cluster_df.to_csv("hotspots.csv", index=False)
    print("[INFO] Saved hotspots.csv")

    # =========================
    # SELECT TOP PER FEATURE
    # =========================
    cluster_df_sdf = (
        cluster_df
        .sort_values(by=["count_points"], ascending=False)
        .groupby("feature", group_keys=False)
        .head(MAX_PER_FEATURE)
        .reset_index(drop=True)
    )

    print(f"[INFO] Hotspots selected for SDF: {len(cluster_df_sdf)}")

    # Create SDF
    hotspot_mol = create_hotspot_mol(cluster_df_sdf)

    writer = Chem.SDWriter("hotspots.sdf")
    writer.write(hotspot_mol)
    writer.close()

    combined = mols + [hotspot_mol]
    writer = Chem.SDWriter("ligands_plus_hotspots.sdf")
    for mol in combined:
        writer.write(mol)
    writer.close()

    print("[DONE] Pharmacophore hotspot analysis complete.")


if __name__ == "__main__":
    main()