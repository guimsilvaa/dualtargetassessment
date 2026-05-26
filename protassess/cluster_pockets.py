import os
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import (linkage, dendrogram)
from sklearn.preprocessing import (StandardScaler)
from sklearn.decomposition import PCA

FEATURE_COLUMNS = [
    "volume",
    "surface_area",
    "surface_volume_ratio",
    "compactness",
    "sphericity",
    "avg_depth",
    "max_depth",
    "avg_hydropathy",
    "hydrophobic_fraction",
    "aromatic_fraction",
    "charged_fraction",
    "n_interface_residues"
]

def cluster_pockets(descriptors_csv, output_dir):
    """
    Perform clustering
    of cavity descriptor profiles.
    """
    logging.info("\nRunning clustering...")

    # Load descriptors
    df = pd.read_csv(
        descriptors_csv
    )
    # Feature matrix
    X = df[
        FEATURE_COLUMNS
    ].values

    # Standardization
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Hierarchical clustering
    linkage_matrix = linkage(
        X_scaled,
        method="ward"
    )
    # Create labels
    labels = [
        f"{target}_{structure}"
        for target, structure in zip(
            df["target"],
            df["structure"]
        )
    ]
    # Plot dendrogram
    plt.figure(figsize=(16, 8))
    dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=8
    )
    plt.title(
        "Hierarchical Clustering Dendrogram"
    )
    plt.ylabel(
        "Ward Distance"
    )
    plt.tight_layout()
    dendrogram_file = os.path.join(
        output_dir,
        "plot_dendrogram.png"
    )
    plt.savefig(
        dendrogram_file,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
    logging.info(
        f"Dendrogram saved: "
        f"{dendrogram_file}"
    )

    # PCA projection
    logging.info(
        "Running PCA projection..."
    )
    pca = PCA(
        n_components=2
    )
    X_pca = pca.fit_transform(
        X_scaled
    )
    df["PC1"] = X_pca[:, 0]
    df["PC2"] = X_pca[:, 1]
    
    # Plot PCA
    plt.figure(figsize=(10, 8))
    targets = sorted(
        df["target"].unique()
    )
    for target in targets:
        subset = df[
            df["target"] == target
        ]
        plt.scatter(
            subset["PC1"],
            subset["PC2"],
            label=target,
            s=70,
            alpha=0.8
        )
    plt.xlabel(
        f"PC1 "
        f"({pca.explained_variance_ratio_[0]*100:.1f}%)"
    )
    plt.ylabel(
        f"PC2 "
        f"({pca.explained_variance_ratio_[1]*100:.1f}%)"
    )
    plt.title(
        "Descriptors PCA"
    )
    plt.legend()
    plt.tight_layout()
    pca_file = os.path.join(
        output_dir,
        "plot_pca.png"
    )
    plt.savefig(
        pca_file,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
    logging.info(
        f"PCA plot saved: "
        f"{pca_file}"
    )

    # PCA loadings
    loading_df = pd.DataFrame(
        pca.components_.T,
        columns=["PC1", "PC2"],
        index=FEATURE_COLUMNS
    )
    loading_file = os.path.join(
        output_dir,
        "plot_pca_feature_loadings.csv"
    )
    loading_df.to_csv(
        loading_file
    )
    logging.info(
        f"PCA loadings saved: "
        f"{loading_file}"
    )

    # Explained variance
    explained_variance = pd.DataFrame({
        "PC": ["PC1", "PC2"],
        "explained_variance_ratio": (
            pca.explained_variance_ratio_
        )
    })
    explained_file = os.path.join(
        output_dir,
        "plot_pca_explained_variance.csv"
    )
    explained_variance.to_csv(
        explained_file,
        index=False
    )
    logging.info(
        f"PCA explained variance saved: "
        f"{explained_file}"
    )
    return df
