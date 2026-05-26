import numpy as np
import pandas as pd
import logging

def compute_residue_composition(residue_list):
    """
    Compute residue class fractions:
    - hydrophobic_fraction
    - aromatic_fraction
    - charged_fraction
    """
    hydrophobic = {
        "ALA", "VAL", "LEU", "ILE",
        "MET", "PRO", "PHE", "TRP"
    }
    aromatic = {
        "PHE", "TYR", "TRP", "HIS"
    }
    charged = {
        "ASP", "GLU", "LYS",
        "ARG", "HIS"
    }
    residue_names = []

    for residue in residue_list:
        # pyKVFinder residue format [resnum, chain, resname]
        if isinstance(residue, (list, tuple)):
            resname = residue[2]
        else:
            tokens = str(residue).split()
            resname = tokens[0]
        residue_names.append(
            resname.upper()
        )
    total = len(residue_names)
    if total == 0:
        return 0.0, 0.0, 0.0
    hydrophobic_fraction = sum(
        r in hydrophobic
        for r in residue_names
    ) / total
    aromatic_fraction = sum(
        r in aromatic
        for r in residue_names
    ) / total
    charged_fraction = sum(
        r in charged
        for r in residue_names
    ) / total
    return (
        hydrophobic_fraction,
        aromatic_fraction,
        charged_fraction
    )

def extract_ligand_proximal_cavity(results, ligand_centroid, debug=False):
    """
    Select cavity closest to ligand centroid
    and compute descriptors.
    """
    cavities_grid = results.cavities
    cavity_ids = list(
        results.volume.keys()
    )
    if len(cavity_ids) == 0:
        raise ValueError(
            "No cavities detected."
        )
    origin = results._vertices[0]
    step = results._step
    cavity_centroids = {}

    # Compute cavity centroids
    for idx, cavity_id in enumerate(cavity_ids):
        voxel_label = idx + 2
        voxel_indices = np.argwhere(
            cavities_grid == voxel_label
        )
        if len(voxel_indices) == 0:
            continue
        xyz_coords = origin + (
            voxel_indices * step
        )
        centroid = xyz_coords.mean(axis=0)
        cavity_centroids[
            cavity_id
        ] = centroid

    # Compute distances to ligand  
    distances = {}
    for cavity_id, centroid in (
        cavity_centroids.items()
    ):
        if np.any(np.isnan(centroid)):
            continue
        dist = np.linalg.norm(
            centroid - ligand_centroid
        )
        distances[cavity_id] = dist
    if len(distances) == 0:
        raise ValueError(
            "No valid cavity distances calculated."
        )

    # Select closest cavity
    selected_cavity = min(
        distances,
        key=distances.get
    )
 
    # Debug logging
    if True:
        logging.info(
            "\nDescriptor-supported cavities:"
        )
        logging.info(cavity_ids)
        logging.info(
            "\nComputed cavity centroids:"
        )
        for cid, cent in (
            cavity_centroids.items()
        ):
            logging.info(
                f"{cid}: {cent}"
            )
        logging.info(
            "\nDistances to ligand:"
        )
        for cid, dist in (
            distances.items()
        ):
            logging.info(
                f"{cid}: {dist:.4f}"
            )
        logging.info(
            f"\nSelected cavity: "
            f"{selected_cavity}"
        )
    
    # Surface descriptors
    surface_area = float(
        results.area[selected_cavity]
    )
    volume = float(
        results.volume[selected_cavity]
    )
    surface_volume_ratio = (
        surface_area / volume
    )
    compactness = (
        volume / surface_area
    )
    sphericity = (
        (np.pi ** (1/3))
        * ((6 * volume) ** (2/3))
    ) / surface_area

    # Residue composition
    (
        hydrophobic_fraction,
        aromatic_fraction,
        charged_fraction
    ) = compute_residue_composition(
        results.residues[
            selected_cavity
        ]
    )

    # Descriptor dictionary
    cavity_data = {
        "cavity_id":
            str(selected_cavity),
        "distance_to_ligand":
            float(
                distances[selected_cavity]
            ),
        "volume":
            volume,
        "surface_area":
            surface_area,
        "surface_volume_ratio":
            surface_volume_ratio,
        "compactness":
            compactness,
        "sphericity":
            sphericity,
        "avg_depth":
            float(
                results.avg_depth[
                    selected_cavity
                ]
            ),
        "max_depth":
            float(
                results.max_depth[
                    selected_cavity
                ]
            ),
        "avg_hydropathy":
            float(
                results.avg_hydropathy[
                    selected_cavity
                ]
            ),
        "hydrophobic_fraction":
            hydrophobic_fraction,
        "aromatic_fraction":
            aromatic_fraction,
        "charged_fraction":
            charged_fraction,
        "n_interface_residues":
            int(
                len(
                    results.residues[
                        selected_cavity
                    ]
                )
            ),
    }
    return (
        cavity_data,
        selected_cavity
    )

def build_descriptor_dataframe(all_results):
    """
    Build final descriptor dataframe.
    """
    df = pd.DataFrame(all_results)
    columns = [
        "target",
        "structure",
        "selected_chain",
        "cavity_id",
        "distance_to_ligand",
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
        "n_interface_residues",
    ]
    return df[columns]