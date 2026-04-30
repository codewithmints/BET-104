import os
import sys
import gzip
import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial
from Bio.PDB import PDBParser
from tqdm import tqdm

# Mapping of amino acid codes to their respective physical size categories
RESIDUE_DIMENSIONS = {
    "G": "Tiny", "A": "Tiny", 
    "V": "Small", "P": "Small", "S": "Small", "T": "Small", "C": "Small",
    "D": "Intermediate", "L": "Intermediate", "I": "Intermediate", "N": "Intermediate",
    "K": "Large", "E": "Large", "M": "Large", "Q": "Large", "H": "Large",
    "R": "Bulky", "F": "Bulky", "Y": "Bulky", "W": "Bulky"
}

def calculate_rotation_angle(vec_a, vec_b, ref_axis, use_degrees=True):
    """
    Computes the signed angle between two vectors around a specific axis.
    """
    # Normalize input vectors
    unit_a = vec_a / np.linalg.norm(vec_a)
    unit_b = vec_b / np.linalg.norm(vec_b)
    unit_axis = ref_axis / np.linalg.norm(ref_axis)

    # Vector geometry calculations
    orthogonal_comp = np.cross(unit_a, unit_b)
    scalar_comp = np.dot(unit_a, unit_b)

    # Determine signed direction
    rad_angle = np.arctan2(np.dot(orthogonal_comp, unit_axis), scalar_comp)

    return np.degrees(rad_angle) if use_degrees else rad_angle

def fetch_structure_data(code, directory):
    """Safely unzips and parses a PDB file."""
    src_path = os.path.join(directory, f"{code}.pdb.gz")
    if not os.path.exists(src_path):
        return None
    
    pdb_engine = PDBParser(QUIET=True)
    try:
        with gzip.open(src_path, "rt") as handle:
            return pdb_engine.get_structure(code, handle)
    except Exception:
        return None

def extract_specific_residue(entity, chain_name, seq_index):
    """Locates a residue by its sequence number within a chain."""
    try:
        target_chain = entity[0][chain_name] # Look in first model
        for item in target_chain:
            if item.id[1] == seq_index:
                return item
    except (KeyError, Exception):
        return None
    return None

def locate_sidechain_centroid(residue_obj):
    """Calculates the center of mass for sidechain atoms only."""
    backbone_atoms = {"N", "CA", "C", "O"}
    atom_positions = [a.get_coord() for a in residue_obj if a.get_name() not in backbone_atoms]
    
    if not atom_positions:
        return None
    return np.mean(atom_positions, axis=0)

def analyze_structural_context(tsv_name, input_path, pdb_path, target_residue_type):
    """
    Parses structural TSVs to compute geometric angles between residues.
    """
    full_path = os.path.join(input_path, tsv_name)
    if os.path.getsize(full_path) == 0:
        return [], None

    try:
        dataset = pd.read_csv(full_path, sep="\t", header=None)
    except Exception:
        return [], None

    if dataset.empty:
        return [], None

    # Identify the PDB and load its 3D coordinates
    current_pdb = dataset.iloc[0, 12]
    protein_model = fetch_structure_data(current_pdb, pdb_path)
    if protein_model is None:
        return [], None

    extracted_results = []
    # Iterate through chunks of 3 (previous, current, next)
    for idx in range(0, len(dataset), 3):
        try:
            prev_row, curr_row, next_row = dataset.iloc[idx], dataset.iloc[idx+1], dataset.iloc[idx+2]
        except IndexError:
            continue

        # Filter for secondary structure and specific AA target
        if curr_row[0] != target_residue_type or curr_row[11] != "HHH":
            continue

        c_id = curr_row[1]
        # Map IDs to actual PDB objects
        r_1 = extract_specific_residue(protein_model, c_id, int(prev_row[2]))
        r_2 = extract_specific_residue(protein_model, c_id, int(curr_row[2]))
        r_3 = extract_specific_residue(protein_model, c_id, int(next_row[2]))

        if not all([r_1, r_2, r_3]):
            continue

        # Extract Alpha-Carbon coordinates
        coords_ca = {
            'p': r_1["CA"].get_coord() if "CA" in r_1 else None,
            'c': r_2["CA"].get_coord() if "CA" in r_2 else None,
            'n': r_3["CA"].get_coord() if "CA" in r_3 else None
        }
        
        # Extract Centroids
        cent_p = locate_sidechain_centroid(r_1)
        cent_c = locate_sidechain_centroid(r_2)

        if any(v is None for v in [coords_ca['p'], coords_ca['c'], coords_ca['n'], cent_p, cent_c]):
            continue

        # Compute geometric vectors
        vector_prev = cent_p - coords_ca['p']
        vector_curr = cent_c - coords_ca['c']
        rotation_axis = coords_ca['c'] - coords_ca['p']

        if np.linalg.norm(rotation_axis) == 0:
            continue

        try:
            measured_angle = calculate_rotation_angle(vector_prev, vector_curr, rotation_axis)
        except Exception:
            continue

        res_type = prev_row[9]
        category = RESIDUE_DIMENSIONS.get(res_type, "Unknown")
        
        if category != "Unknown":
            extracted_results.append([current_pdb, res_type, category, measured_angle])

    return extracted_results, (current_pdb if extracted_results else None)

def run_main_workflow():
    # CLI Argument unpacking
    src_dir = sys.argv[1]
    save_path = sys.argv[2]
    query_aa = sys.argv[3]
    pdb_repo = sys.argv[4] if len(sys.argv) > 4 else "pdbs"
    thread_count = int(sys.argv[5]) if len(sys.argv) > 5 else 12

    context_files = sorted([f for f in os.listdir(src_dir) if f.endswith(".tsv")])

    # Curry the processing function with static paths
    task_runner = partial(
        analyze_structural_context,
        input_path=src_dir,
        pdb_path=pdb_repo,
        target_residue_type=query_aa
    )

    collected_data = []
    successful_pdbs = set()

    # Parallel execution block
    with mp.Pool(thread_count) as executor:
        execution_iterator = executor.imap_unordered(task_runner, context_files, chunksize=16)
        
        for records, pdb_name in tqdm(execution_iterator, total=len(context_files), desc="Analyzing Proteins"):
            if records:
                collected_data.extend(records)
            if pdb_name:
                successful_pdbs.add(str(pdb_name))

    # Output assembly
    final_df = pd.DataFrame(collected_data, columns=["pdb", "left_aa", "size_class", "angle"])
    base_folder = os.path.dirname(save_path) or "."
    os.makedirs(base_folder, exist_ok=True)
    
    final_df.to_csv(save_path, sep="\t", index=False)

    # Save list of processed PDBs
    pdb_log = os.path.join(base_folder, "valid_pdbs.txt")
    with open(pdb_log, "w") as log_file:
        for item in sorted(successful_pdbs):
            log_file.write(f"{item}\n")

    print(f"Angle data written to: {save_path}")
    print(f"Log of valid PDBs: {pdb_log}")
    print(f"Entries processed: {len(final_df)}")

if __name__ == "__main__":
    run_main_workflow()
