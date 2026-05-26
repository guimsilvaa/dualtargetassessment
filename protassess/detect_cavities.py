import os
import logging
import pickle
import warnings
import pyKVFinder
import numpy as np
from tqdm import tqdm
from utils import (clean_pdb_structure, get_ligand_centroid)
from extract_descriptors import (extract_ligand_proximal_cavity)

warnings.filterwarnings("ignore", category=UserWarning)
DEBUG = True

def export_cavity_pdb(results, cavity_id, output_pdb):
    """
    Export selected cavity voxels
    as pseudoatom PDB file.
    """
    cavities_grid = results.cavities
    origin = results._vertices[0]
    step = results._step
    
    # Map descriptor cavity labels to voxel grid labels
    descriptor_ids = list(
        results.volume.keys()
    )
    if cavity_id not in descriptor_ids:
        raise ValueError(
            f"Cavity ID {cavity_id} "
            f"not found in descriptor list."
        )
    voxel_label = (
        descriptor_ids.index(cavity_id) + 2
    )

    # Extract voxel coordinates
    voxel_indices = np.argwhere(
        cavities_grid == voxel_label
    )
    if len(voxel_indices) == 0:
        raise ValueError(
            f"No voxels found for cavity "
            f"{cavity_id}"
        )

    # Convert to Cartesian coordinates
    xyz_coords = origin + (
        voxel_indices * step
    )

    # Write pseudoatom PDB
    with open(output_pdb, "w") as f:
        f.write(
            f"REMARK Cavity ID: "
            f"{cavity_id}\n"
        )
        f.write(
            f"REMARK Number of voxels: "
            f"{len(xyz_coords)}\n"
        )
        for idx, coord in enumerate(
            xyz_coords
        ):
            x, y, z = coord
            line = (
                f"HETATM"
                f"{idx+1:5d} "
                f" C   "
                f"CAV A   1    "
                f"{x:8.3f}"
                f"{y:8.3f}"
                f"{z:8.3f}"
                f"  1.00  1.00           C\n"
            )
            f.write(line)
        f.write("END\n")
    logging.info(
        f"Cavity exported: "
        f"{output_pdb}"
    )

def process_target_folder(target_folder, output_dir ):
    """
    Process one target ensemble.
    """
    logging.info(
        f"\nProcessing target: "
        f"{target_folder}"
    )
    # Output folders
    cleaned_dir = os.path.join(
        output_dir,
        f"{target_folder}_cleaned_pdbs"
    )
    visualization_dir = os.path.join(
        output_dir,
        f"{target_folder}_cavity_visualizations"
    )
    os.makedirs(cleaned_dir, exist_ok=True)
    os.makedirs(
        visualization_dir,
        exist_ok=True
    )
    # Locate ligand SDF
    sdf_files = [
        f for f in os.listdir(
            target_folder
        )
        if f.endswith(".sdf")
    ]
    if len(sdf_files) != 1:
        raise ValueError(
            f"{target_folder} must contain "
            f"exactly ONE reference ligand SDF."
        )
    ligand_sdf = os.path.join(
        target_folder,
        sdf_files[0]
    )
    # Ligand centroid
    ligand_centroid = (
        get_ligand_centroid(
            ligand_sdf
        )
    )
    logging.info(
        f"Ligand centroid: "
        f"{ligand_centroid}"
    )
    # Locate PDB files
    pdb_files = sorted([
        f for f in os.listdir(
            target_folder
        )
        if f.endswith(".pdb")
    ])
    if len(pdb_files) == 0:
        raise ValueError(
            f"No PDB files found in "
            f"{target_folder}"
        )
    all_results = []

    # Main loop
    for pdb_file in tqdm(
        pdb_files,
        desc=target_folder
    ):
        structure_name = pdb_file
        input_pdb = os.path.join(
            target_folder,
            pdb_file
        )
        cleaned_pdb = os.path.join(
            cleaned_dir,
            pdb_file.replace(
                ".pdb",
                "_cleaned.pdb"
            )
        )
        try:
            # Clean structure
            cleaned_pdb, selected_chain = (
                clean_pdb_structure(
                    input_pdb=input_pdb,
                    output_pdb=cleaned_pdb,
                    ligand_sdf=ligand_sdf
                )
            )
            # Run pyKVFinder
            results = pyKVFinder.run_workflow(
                cleaned_pdb,
                include_depth=True,
                include_hydropathy=True,
                ignore_backbone=True,
            )
            # Extract descriptors
            cavity_data, selected_cavity = (
                extract_ligand_proximal_cavity(
                    results,
                    ligand_centroid,
                    debug=DEBUG
                )
            )
            # Export cavity visualization
            cavity_pdb = os.path.join(
                visualization_dir,
                pdb_file.replace(
                    ".pdb",
                    "_cavity.pdb"
                )
            )
            export_cavity_pdb(
                results,
                selected_cavity,
                cavity_pdb
            )
            # Metadata
            cavity_data["target"] = (
                target_folder
            )
            cavity_data["structure"] = (
                structure_name
            )
            cavity_data["selected_chain"] = (
                selected_chain
            )
            cavity_data["cleaned_pdb"] = (
                cleaned_pdb
            )
            cavity_data["cavity_pdb"] = (
                cavity_pdb
            )
            all_results.append(
                cavity_data
            )
            logging.info(
                f"SUCCESS: {structure_name}"
            )
        except Exception as e:
            logging.exception(
                f"ERROR processing "
                f"{structure_name}: {e}"
            )
    return all_results