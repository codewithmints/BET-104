import sys
import os
import re
import pandas as pd

# Mapping for 3-letter to 1-letter amino acid codes
RESIDUE_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

def build_structural_context():
    # Input argument handling
    input_ss = sys.argv[1]
    save_path = sys.argv[2]
    aa_query = sys.argv[3]

    # Extracting filename without extension for tracking
    source_tag = os.path.basename(input_ss).split(os.extsep)[0]

    parsed_data = []

    # Stream the file and filter for Assignment (ASG) lines
    with open(input_ss, 'r') as stream:
        for entry in stream:
            if not entry.startswith("ASG"):
                continue
            
            elements = entry.split()
            if len(elements) < 10:
                continue

            # Identify the sequence number using regex
            num_match = re.search(r"\d+", elements[3])
            if not num_match:
                continue

            try:
                # Store structural attributes in a flat list for efficiency
                parsed_data.append([
                    elements[1],           # 3-letter code
                    elements[2],           # Chain ID
                    int(num_match.group()),# Sequence number
                    int(elements[4]),      # Index
                    elements[5],           # Secondary structure code
                    elements[6],           # Secondary structure name
                    float(elements[7]),    # Phi angle
                    float(elements[8]),    # Psi angle
                    float(elements[9])     # Solvent accessibility area
                ])
            except (ValueError, IndexError):
                continue

    if not parsed_data:
        return

    # Convert to DataFrame for sliding window analysis
    table = pd.DataFrame(parsed_data, columns=[
        "res3", "chain", "resnum", "idx", "ss_code", "ss_name", "phi", "psi", "area"
    ])

    triplet_results = []

    # Iterate through the sequence to find the target amino acid and its neighbors
    for idx in range(1, len(table) - 1):
        mid = table.iloc[idx]
        
        # Only process if the center matches our target residue
        if mid["res3"] != aa_query:
            continue

        pre = table.iloc[idx - 1]
        aft = table.iloc[idx + 1]

        # Generate the single-letter triplet string (e.g., "HHH")
        aa_sequence = "".join([
            RESIDUE_MAP.get(pre["res3"], "X"),
            RESIDUE_MAP.get(mid["res3"], "X"),
            RESIDUE_MAP.get(aft["res3"], "X")
        ])

        # Combine secondary structure codes
        ss_sequence = pre["ss_code"] + mid["ss_code"] + aft["ss_code"]

        # Formulate a location identifier string
        loc_map = f"{mid['chain']}:{pre['resnum']},{mid['chain']}:{mid['resnum']},{mid['chain']}:{aft['resnum']}"

        # Flatten the neighbors into the final output format
        for member in [pre, mid, aft]:
            triplet_results.append([
                member["res3"],
                member["chain"],
                member["resnum"],
                member["idx"],
                member["ss_code"],
                member["ss_name"],
                member["phi"],
                member["psi"],
                member["area"],
                RESIDUE_MAP.get(member["res3"], "X"),
                aa_sequence,
                ss_sequence,
                source_tag,
                loc_map
            ])

    # Ensure output directory exists before saving
    target_dir = os.path.dirname(save_path)
    if target_dir and not os.path.exists(target_dir):
        os.makedirs(target_dir, exist_ok=True)

    # Save results as a tab-separated file without headers
    final_output = pd.DataFrame(triplet_results)
    final_output.to_csv(save_path, sep="\t", index=False, header=False)

if __name__ == "__main__":
    build_structural_context()
