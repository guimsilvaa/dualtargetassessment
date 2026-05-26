import logging
import numpy as np
from Bio.PDB import (PDBParser, PDBIO, Select, is_aa)
from rdkit import Chem

def get_ligand_centroid(sdf_file):
    """
    Compute centroid coordinates
    from reference ligand SDF.
    """
    mol = Chem.SDMolSupplier(
        sdf_file,
        removeHs=False
    )[0]
    conf = mol.GetConformer()
    coords = np.array([
        list(
            conf.GetAtomPosition(i)
        )
        for i in range(
            mol.GetNumAtoms()
        )
    ])
    centroid = coords.mean(axis=0)
    return centroid

def detect_closest_chain(structure, ligand_centroid):
    """
    Detect protein chain closest
    to ligand centroid.
    """
    chain_distances = {}
    for model in structure:
        for chain in model:
            protein_coords = []
            for residue in chain:
                if not is_aa(residue):
                    continue
                for atom in residue:
                    protein_coords.append(
                        atom.coord
                    )
            if len(protein_coords) == 0:
                continue
            protein_coords = np.array(
                protein_coords
            )
            distances = np.linalg.norm(
                protein_coords
                - ligand_centroid,
                axis=1
            )
            min_distance = np.min(
                distances
            )
            chain_distances[
                chain.id
            ] = min_distance
    if len(chain_distances) == 0:
        raise ValueError(
            "No protein chains detected."
        )

    # Select closest chain
    selected_chain = min(
        chain_distances,
        key=chain_distances.get
    )

    # Logging
    logging.info(
        "\nDetected chain distances to ligand:"
    )
    for cid, dist in (
        chain_distances.items()
    ):
        logging.info(
            f"Chain {cid}: "
            f"{dist:.2f} Å"
        )
    logging.info(
        f"Selected chain: "
        f"{selected_chain}"
    )
    return selected_chain

# CUSTOM BIOPYTHON SELECTOR
class ProteinChainSelect(Select):
    def __init__(
        self,
        selected_chain
    ):
        self.selected_chain = (
            selected_chain
        )

    def accept_chain(
        self,
        chain
    ):
        return (
            chain.id
            == self.selected_chain
        )

    def accept_residue(
        self,
        residue
    ):
        # Keep only amino acids
        return is_aa(residue)

    def accept_atom(
        self,
        atom
    ):
        # Remove hydrogens
        if atom.element == "H":
            return False

        # Remove alternate conformers
        altloc = atom.get_altloc()
        if altloc not in (
            " ",
            "",
            "A"
        ):
            return False
        return True


def clean_pdb_structure(input_pdb, output_pdb, ligand_sdf):
    """
    Cleaning procedure:
    - detect closest protein chain
    - keep only amino acids
    - remove waters
    - remove ligands
    - remove ions/metals
    - remove alternate conformers
    - keep single monomer chain
    """
    logging.info(
        f"\nProcessing structure: "
        f"{input_pdb}"
    )

    # Parse structure
    parser = PDBParser(
        QUIET=True
    )
    structure = parser.get_structure(
        "protein",
        input_pdb
    )

    # Ligand centroid
    ligand_centroid = (
        get_ligand_centroid(
            ligand_sdf
        )
    )

    # Detect closest chain
    selected_chain = (
        detect_closest_chain(
            structure,
            ligand_centroid
        )
    )

    # Export cleaned structure
    io = PDBIO()
    io.set_structure(
        structure
    )
    io.save(
        output_pdb,
        ProteinChainSelect(
            selected_chain
        )
    )
    logging.info(
        f"Cleaned structure saved: "
        f"{output_pdb}"
    )
    return (
        output_pdb,
        selected_chain
    )