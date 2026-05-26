import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import (linkage, leaves_list)

DESCRIPTORS = [
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

def compute_cohens_d(x1, x2):
    mean1 = np.mean(x1)
    mean2 = np.mean(x2)

    std1 = np.std(x1, ddof=1)
    std2 = np.std(x2, ddof=1)

    n1 = len(x1)
    n2 = len(x2)

    pooled_std = np.sqrt(
        (
            ((n1 - 1) * std1**2)
            +
            ((n2 - 1) * std2**2)
        )
        /
        (n1 + n2 - 2)
    )
    d = (mean1 - mean2) / pooled_std
    return d

def analyze_ensembles(descriptors_csv, output_dir):
    logging.info(
        "\nRunning ensemble analysis..."
    )

    # Load descriptor table
    df = pd.read_csv(
        descriptors_csv
    )
    targets = sorted(
        df["target"].unique()
    )
    if len(targets) != 2:
        raise ValueError(
            "Current implementation "
            "supports exactly 2 targets."
        )
    target1 = targets[0]
    target2 = targets[1]
    df1 = df[
        df["target"] == target1
    ]
    df2 = df[
        df["target"] == target2
    ]

    # Output folders
    plots_dir = os.path.join(
        output_dir,
        "plots"
    )
    os.makedirs(
        plots_dir,
        exist_ok=True
    )
    summary_results = []
    importance_results = []

    # VIOLIN PLOTS + STATS
    for descriptor in DESCRIPTORS:
        logging.info(
            f"Analyzing descriptor: "
            f"{descriptor}"
        )
        values1 = df1[
            descriptor
        ].dropna()
        values2 = df2[
            descriptor
        ].dropna()

        # Statistics
        t_stat, p_value = ttest_ind(
            values1,
            values2,
            equal_var=False
        )
        cohens_d = compute_cohens_d(
            values1,
            values2
        )
        summary_results.append({
            "descriptor":
                descriptor,
            "target1_mean":
                np.mean(values1),
            "target1_std":
                np.std(values1),
            "target2_mean":
                np.mean(values2),
            "target2_std":
                np.std(values2),
            "pvalue":
                p_value,
            "cohens_d":
                cohens_d
        })
        importance_results.append({
            "descriptor":
                descriptor,
            "cohens_d":
                abs(cohens_d),
            "signed_cohens_d":
                cohens_d,
            "pvalue":
                p_value
        })

        # Violin plot
        plt.figure(figsize=(7, 6))
        violin_data = [
            values1,
            values2
        ]
        parts = plt.violinplot(
            violin_data,
            showmeans=True,
            showextrema=True
        )
        colors = [ 
                "#1f77b4", # blue 
                "#ff7f0e" # orange 
        ] 
        for patch, color in zip( 
                parts["bodies"], 
                colors 
        ): 
                patch.set_facecolor(color) 
                patch.set_edgecolor("black") 
                patch.set_alpha(0.7)
        plt.xticks(
            [1, 2],
            [target1, target2]
        )
        plt.ylabel(
            descriptor
        )
        plt.title(
            f"{descriptor} distribution"
        )
        plt.tight_layout()
        plot_file = os.path.join(
            plots_dir,
            f"{descriptor}_violin.png"
        )
        plt.savefig(
            plot_file,
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()

    # SUMMARY CSV
    summary_df = pd.DataFrame(
        summary_results
    )
    summary_csv = os.path.join(
        output_dir,
        "results_summary.csv"
    )
    summary_df.to_csv(
        summary_csv,
        index=False
    )
    logging.info(
        f"Summary CSV saved: "
        f"{summary_csv}"
    )

    # DESCRIPTOR IMPORTANCE
    importance_df = pd.DataFrame(
        importance_results
    )
    importance_df = importance_df.sort_values(
        by="cohens_d",
        ascending=False
    )
    importance_csv = os.path.join(
        output_dir,
        "descriptor_importance.csv"
    )
    importance_df.to_csv(
        importance_csv,
        index=False
    )
    logging.info(
        f"Descriptor importance saved: "
        f"{importance_csv}"
    )

    # Importance plot
    plt.figure(figsize=(10, 6))
    plt.barh(
        importance_df["descriptor"],
        importance_df["cohens_d"]
    )
    plt.xlabel(
        "Absolute Cohen's d"
    )
    plt.title(
        "Descriptor Importance Ranking"
    )
    plt.gca().invert_yaxis()
    plt.tight_layout()
    importance_plot = os.path.join(
        output_dir,
        "plot_descriptor_importance.png"
    )
    plt.savefig(
        importance_plot,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # HEATMAP
    logging.info(
        "Generating descriptor heatmap..."
    )
    heatmap_df = df.copy()
    X = heatmap_df[
        DESCRIPTORS
    ].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    scaled_df = pd.DataFrame(
        X_scaled,
        columns=DESCRIPTORS
    )
    scaled_df["target"] = (
        heatmap_df["target"].values
    )
    scaled_df["structure"] = (
        heatmap_df["structure"].values
    )
    
    # Cluster rows
    linkage_matrix = linkage(
        X_scaled,
        method="ward"
    )
    ordered_indices = leaves_list(
        linkage_matrix
    )
    scaled_matrix = scaled_df[
        DESCRIPTORS
    ].iloc[
        ordered_indices
    ]
    labels = [
        f"{t}_{s}"
        for t, s in zip(
            scaled_df["target"].iloc[
                ordered_indices
            ],
            scaled_df["structure"].iloc[
                ordered_indices
            ]
        )
    ]
    # Plot heatmap
    plt.figure(figsize=(14, 10))
    plt.imshow(
        scaled_matrix,
        aspect="auto",
        cmap="coolwarm"
    )
    plt.colorbar(
        label="Z-score"
    )
    plt.yticks(
        np.arange(len(labels)),
        labels,
        fontsize=6
    )
    plt.xticks(
        np.arange(len(DESCRIPTORS)),
        DESCRIPTORS,
        rotation=90
    )
    plt.title(
        "Descriptor Heatmap"
    )
    plt.tight_layout()
    heatmap_file = os.path.join(
        output_dir,
        "plot_heatmap.png"
    )
    plt.savefig(
        heatmap_file,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # RADAR PLOT
    logging.info(
        "Generating radar plot..."
    )
    radar_df = pd.DataFrame(
        X_scaled,
        columns=DESCRIPTORS
    )
    radar_df["target"] = (
        df["target"].values
    )
    mean_profiles = radar_df.groupby(
        "target"
    ).mean()
    categories = DESCRIPTORS
    N = len(categories)
    angles = np.linspace(
        0,
        2 * np.pi,
        N,
        endpoint=False
    ).tolist()
    angles += angles[:1]
    fig, ax = plt.subplots(
        figsize=(8, 8),
        subplot_kw=dict(
            polar=True
        )
    )
    for target in mean_profiles.index:
        values = mean_profiles.loc[
            target
        ].tolist()
        values += values[:1]
        ax.plot(
            angles,
            values,
            linewidth=2,
            label=target
        )
        ax.fill(
            angles,
            values,
            alpha=0.2
        )
    ax.set_xticks(
        angles[:-1]
    )
    ax.set_xticklabels(
        categories,
        fontsize=8
    )
    plt.title(
        "Mean Descriptor Profile"
    )
    plt.legend(
        loc="upper left",
        bbox_to_anchor=(1.15, 1.1)
    )
    radar_file = os.path.join(
        output_dir,
        "plot_radar_plot.png"
    )
    plt.savefig(
        radar_file,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
    logging.info(
        "Ensemble analysis complete."
    )
    return summary_df