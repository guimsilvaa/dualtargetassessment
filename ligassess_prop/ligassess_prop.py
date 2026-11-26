#!/usr/bin/env python3
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy.stats import mannwhitneyu, ks_2samp
# ===============================
PROT1_NAME = "AChE"
PROT2_NAME = "GSK3b"
# ===============================
#STANDARDIZATION
def standardize_mol(mol):
    try:
        mol = rdMolStandardize.FragmentParent(mol)
        normalizer = rdMolStandardize.Normalizer()
        mol = normalizer.normalize(mol)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None
# ===============================
# LOAD SDF
def load_sdf(file, standardize=True):
    suppl = Chem.SDMolSupplier(file, removeHs=False)
    mols = []

    for mol in suppl:
        if mol is None:
            continue

        if standardize:
            mol = standardize_mol(mol)
            if mol is None:
                continue

        mols.append(mol)

    print(f"[INFO] Loaded {len(mols)} valid molecules from {file}")
    return mols
# ===============================
#  PROPERTY CALCULATION
def compute_properties(mol):
    try:
        ring_info = mol.GetRingInfo()
        ring_sizes = [len(r) for r in ring_info.AtomRings()]
        largest_ring = max(ring_sizes) if ring_sizes else 0

        heavy_atoms = Lipinski.HeavyAtomCount(mol)
        rot_bonds = Lipinski.NumRotatableBonds(mol)

        return {
            "MW": Descriptors.MolWt(mol),
            "logP": Descriptors.MolLogP(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "RotBonds": rot_bonds,
            "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
            "HeavyAtoms": heavy_atoms,
            "RingCount": Lipinski.RingCount(mol),
            "AromaticRings": Lipinski.NumAromaticRings(mol),
            "QED": QED.qed(mol),
            "FormalCharge": Chem.GetFormalCharge(mol),
            "HeteroAtoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
            "LargestRingSize": largest_ring,
            "FlexibilityIndex": rot_bonds / heavy_atoms if heavy_atoms > 0 else 0,
        }

    except Exception:
        return None
# ===============================
#  PROCESS DATASET
def process_dataset(file, protein_label, dataset_label):
    mols = load_sdf(file)
    data = []

    for idx, mol in enumerate(mols):
        props = compute_properties(mol)
        if props is None:
            continue

        entry = {
            "molecule_id": f"{protein_label}_{dataset_label}_{idx+1}",
            "protein": protein_label,
            "dataset": dataset_label,
        }

        entry.update(props)
        data.append(entry)

    return pd.DataFrame(data)
# ===============================
# PLOTTING
def plot_property(df, prop, outdir):

    plt.figure(figsize=(7, 5))

    # Define color palette
    palette = {
        f"{PROT1_NAME}_PDB": "#9ecae1",     # light blue
        f"{PROT1_NAME}_ChEMBL": "#1f77b4",  # darker blue
        f"{PROT2_NAME}_PDB": "#fdd0a2",     # light orange
        f"{PROT2_NAME}_ChEMBL": "#e6550d",  # darker orange
    }

    groups = df["group"].unique()

    for g in groups:
        subset = df[df["group"] == g]

        sns.kdeplot(
            data=subset,
            x=prop,
            label=g,
            color=palette.get(g, None),
            fill=True,
            alpha=0.3,
            linewidth=2
        )

    plt.title(f"Distribution of {prop}")
    plt.xlabel(prop)
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()

    outfile = os.path.join(outdir, f"{prop}_density.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
# ===============================
# CLIFF'S DELTA
def cliffs_delta(x, y):
    x = np.array(x)
    y = np.array(y)

    n = len(x)
    m = len(y)

    greater = sum(np.sum(xi > y) for xi in x)
    less = sum(np.sum(xi < y) for xi in x)

    return (greater - less) / (n * m)
# ===============================
def plot_heatmap(stats_df, outdir):

    # Keep only cross-protein comparisons
    mask = stats_df["Comparison"].str.contains("_PDB vs .*_PDB") | \
           stats_df["Comparison"].str.contains("_ChEMBL vs .*_ChEMBL")

    df = stats_df[mask].copy()

    # Simplify labels
    def simplify_label(comp):
        if "PDB" in comp:
            return "PDB"
        elif "ChEMBL" in comp:
            return "ChEMBL"
        return comp

    df["Comparison"] = df["Comparison"].apply(simplify_label)

    # Pivot tables
    pivot = df.pivot(index="Property", columns="Comparison", values="Cliffs_delta")
    sig_mask = df.pivot(index="Property", columns="Comparison", values="Significant")

    # Sort by strongest effect
    pivot["max_abs"] = pivot.abs().max(axis=1)
    pivot = pivot.sort_values("max_abs", ascending=False)
    pivot = pivot.drop(columns="max_abs")

    sig_mask = sig_mask.loc[pivot.index]

    # Plot
    plt.figure(figsize=(5, max(6, len(pivot) * 0.4)))

    ax = sns.heatmap(
        pivot,
        cmap="PiYG",          # <-- NEW COLOR SCHEME
        center=0,
        annot=False,          # we'll annotate manually
        linewidths=0.5,
        cbar_kws={"label": "Cliff's delta"}
    )

    # Manual annotations with bold for significant
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):

            value = pivot.iloc[i, j]
            is_sig = sig_mask.iloc[i, j]

            if pd.isna(value):
                continue

            ax.text(
                j + 0.5,
                i + 0.5,
                f"{value:.2f}",
                ha='center',
                va='center',
                fontweight='bold' if is_sig else 'normal',
                color='black'
            )

    plt.title(f"Property differences - {PROT1_NAME} vs {PROT2_NAME}")
    plt.ylabel("Property")
    plt.xlabel("")

    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()

    outfile = os.path.join(outdir, "property_heatmap.png")
    plt.savefig(outfile, dpi=300)
    plt.close()

    print("[INFO] Saved heatmap: property_heatmap.png")
# ===============================
# STATISTICS
def compare_distributions(df, prop, g1, g2):
    x = df[df["group"] == g1][prop].dropna()
    y = df[df["group"] == g2][prop].dropna()
    if len(x) < 5 or len(y) < 5:
        return None
    try:
        mwu_p = mannwhitneyu(x, y, alternative='two-sided').pvalue
        ks_p = ks_2samp(x, y).pvalue
        effect = cliffs_delta(x.values, y.values)
    except Exception:
        return None
    return mwu_p, ks_p, effect
def classify_effect(delta):
    ad = abs(delta)
    if ad < 0.147:
        return "negligible"
    elif ad < 0.33:
        return "small"
    elif ad < 0.474:
        return "medium"
    else:
        return "large"
# ===============================
def main():

    if len(sys.argv) != 5:
        print("\nUsage:")
        print("python ligassess_prop.py pdb_prot1.sdf chembl_prot1.sdf pdb_prot2.sdf chembl_prot2.sdf\n")
        sys.exit(1)

    pdb1, chembl1, pdb2, chembl2 = sys.argv[1:]

    print("\n[INFO] Processing datasets...\n")

    df1 = process_dataset(pdb1, PROT1_NAME, "PDB")
    df2 = process_dataset(chembl1, PROT1_NAME, "ChEMBL")
    df3 = process_dataset(pdb2, PROT2_NAME, "PDB")
    df4 = process_dataset(chembl2, PROT2_NAME, "ChEMBL")
    
    df = pd.concat([df1, df2, df3, df4], ignore_index=True)

    df["group"] = df["protein"] + "_" + df["dataset"]

    # ===============================
    # SAVE RAW DATA
    # ===============================
    df.to_csv("ligand_properties.csv", index=False)
    print("[INFO] Saved ligand_properties.csv")

    # ===============================
    # PLOTS
    # ===============================
    outdir = "plots"
    os.makedirs(outdir, exist_ok=True)

    properties = [
        "MW", "logP", "TPSA", "HBD", "HBA",
        "RotBonds", "FractionCSP3", "HeavyAtoms",
        "RingCount", "AromaticRings", "QED",
        "FormalCharge", "HeteroAtoms",
        "LargestRingSize", "FlexibilityIndex"
    ]
         
    # ===============================
    # STATISTICAL ANALYSIS
    # ===============================
    print("\n[INFO] Running statistical comparisons...\n")

    comparisons = [
        (f"{PROT1_NAME}_PDB", f"{PROT2_NAME}_PDB"),
        (f"{PROT1_NAME}_ChEMBL", f"{PROT2_NAME}_ChEMBL"),
        (f"{PROT1_NAME}_PDB", f"{PROT1_NAME}_ChEMBL"),
        (f"{PROT2_NAME}_PDB", f"{PROT2_NAME}_ChEMBL"),
    ]
    results = []

    for prop in properties:
        for g1, g2 in comparisons:

            res = compare_distributions(df, prop, g1, g2)
            if res is None:
                continue

            mwu_p, ks_p, delta = res
            effect_label = classify_effect(delta)

            significant = (
                (mwu_p < 0.01 or ks_p < 0.01) and abs(delta) >= 0.33
            )

            results.append({
                "Property": prop,
                "Comparison": f"{g1} vs {g2}",
                "MWU_p": mwu_p,
                "KS_p": ks_p,
                "Cliffs_delta": delta,
                "Effect_size": effect_label,
                "Significant": significant
            })

    stats_df = pd.DataFrame(results)
    stats_df.to_csv("ligand_property_stats.csv", index=False)

    # ===============================
    # FILTER SIGNIFICANT RESULTS
    # ===============================
    sig_df = stats_df[stats_df["Significant"] == True]
    sig_df.to_csv("ligand_property_significant.csv", index=False)
    print("[INFO] Saved ligand_property_stats.csv")
    print("[INFO] Saved ligand_property_significant.csv")
    
    plot_heatmap(stats_df, outdir)    
    print("\n[INFO] Generating plots (only significant properties)...\n")
    significant_props = stats_df[stats_df["Significant"] == True]["Property"].unique()
    if len(significant_props) == 0:
        print("[WARNING] No significant properties found. No plots generated.")
    else:
        for prop in significant_props:
            plot_property(df, prop, outdir)
            print(f"[INFO] Plot saved: {prop}_density.png")

    print("\n[DONE] Analysis complete.\n")

if __name__ == "__main__":
    main()