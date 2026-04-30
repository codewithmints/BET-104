"""
Optimized parallel execution engine for structural PDB processing.
Orchestrates file decompression, STRIDE secondary structure assignment, 
and triplet context extraction across multiple CPU cores.
"""

import gzip
import multiprocessing as mp
import os
import re
import subprocess
import sys
import tempfile
import yaml
from pathlib import Path
from tqdm import tqdm
import pandas as pd

# Path configuration
BASE_PATH = Path(__file__).resolve().parent.parent
PDB_STORAGE = BASE_PATH / "pdbs"
STRIDE_CACHE = BASE_PATH / "stride_out"
RESULTS_DIR = BASE_PATH / "contexts"

# Load global configurations
with open(BASE_PATH / "config.yaml") as config_handle:
    SETTINGS = yaml.safe_load(config_handle)

TARGET_RESIDUE = SETTINGS["target_aa"]
STRIDE_EXE = SETTINGS.get("stride_path", "/usr/local/bin/stride")

# Residue nomenclature mapping
MAP_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

def parse_stride_output(filepath):
    """Parses STRIDE output to extract residue-level structural data."""
    structural_records = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith("ASG"):
                continue
            cols = line.split()
            if len(cols) < 10:
                continue
            
            # Numeric extraction for residue sequence
            match = re.match(r"\d+", cols[3])
            if not match:
                continue
            
            try:
                structural_records.append({
                    "res3": cols[1],
                    "chain": cols[2],
                    "resnum": int(match.group()),
                    "idx": int(cols[4]),
                    "ss_code": cols[5],
                    "ss_name": cols[6],
                    "phi": float(cols[7]),
                    "psi": float(cols[8]),
                    "area": float(cols[9]),
                })
            except (ValueError, IndexError):
                continue
    return structural_records

def write_context_triplets(data, query_aa, protein_id, destination):
    """Identifies and saves triplet contexts centered on the query residue."""
    triplet_rows = []
    count = len(data)
    
    for i in range(1, count - 1):
        mid = data[i]
        if mid["res3"] != query_aa:
            continue
            
        pre = data[i - 1]
        aft = data[i + 1]

        # Concatenate single-letter codes and SS codes
        seq_triplet = (
            MAP_3TO1.get(pre["res3"], "X") +
            MAP_3TO1.get(mid["res3"], "X") +
            MAP_3TO1.get(aft["res3"], "X")
        )
        ss_triplet = pre["ss_code"] + mid["ss_code"] + aft["ss_code"]
        
        # Define spatial coordinates string
        coordinate_id = (
            f"{mid['chain']}:{pre['resnum']},"
            f"{mid['chain']}:{mid['resnum']},"
            f"{mid['chain']}:{aft['resnum']}"
        )

        for entry in (pre, mid, aft):
            triplet_rows.append([
                entry["res3"], entry["chain"], entry["resnum"], entry["idx"],
                entry["ss_code"], entry["ss_name"], entry["phi"], entry["psi"], entry["area"],
                MAP_3TO1.get(entry["res3"], "X"), seq_triplet, ss_triplet, protein_id, coordinate_id,
            ])

    if triplet_rows:
        pd.DataFrame(triplet_rows).to_csv(destination, sep="\t", index=False, header=False)

def execute_single_pdb(archive_path):
    """Decompresses, analyzes, and extracts data for a single PDB archive."""
    pdb_name = os.path.basename(archive_path).replace(".pdb.gz", "")
    stride_file = STRIDE_CACHE / f"{pdb_name}.ss.out"
    output_tsv = RESULTS_DIR / f"context_for_{TARGET_RESIDUE}_in_{pdb_name}.tsv"

    if output_tsv.exists():
        return ("already_exists", pdb_name)

    temp_path = None
    try:
        # Create a temporary file for STRIDE to read
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            temp_path = tmp.name
            with gzip.open(archive_path, "rb") as gz_in:
                # Optimized block reading
                for buffer in iter(lambda: gz_in.read(65536), b""):
                    tmp.write(buffer)

        # Run STRIDE external process
        with open(stride_file, "w") as out_log:
            try:
                subprocess.run(
                    [STRIDE_EXE, temp_path],
                    stdout=out_log,
                    stderr=subprocess.DEVNULL,
                    check=False,
                    timeout=60,
                )
            except subprocess.TimeoutExpired:
                return ("timeout", pdb_name)

        # Process the resulting structural data
        extracted_data = parse_stride_output(stride_file)
        write_context_triplets(extracted_data, TARGET_RESIDUE, pdb_name, output_tsv)
        return ("success", pdb_name)

    except Exception as e:
        return (f"fail:{type(e).__name__}", pdb_name)
    finally:
        if temp_path and os.path.exists(temp_path):
            try:
                os.unlink(temp_path)
            except OSError:
                pass

def run_orchestrator():
    """Main entry point for batch processing."""
    concurrency = int(sys.argv[1]) if len(sys.argv) > 1 else 12

    # Ensure necessary folders exist
    STRIDE_CACHE.mkdir(exist_ok=True)
    RESULTS_DIR.mkdir(exist_ok=True)

    target_files = sorted(str(p) for p in PDB_STORAGE.glob("*.pdb.gz"))
    
    print(f"--- PDB Pipeline Driver ---")
    print(f"File Count: {len(target_files)}")
    print(f"Core Count: {concurrency}")
    print(f"Residue:    {TARGET_RESIDUE}")

    stats = {}
    with mp.Pool(concurrency) as pool:
        # Execute pool with progress bar
        job_iterator = pool.imap_unordered(execute_single_pdb, target_files, chunksize=8)
        for outcome, _ in tqdm(job_iterator, total=len(target_files), desc="Running Pipeline"):
            stats[outcome] = stats.get(outcome, 0) + 1

    print("\nExecution Summary:")
    for status, count in sorted(stats.items()):
        print(f"  {status}: {count}")

if __name__ == "__main__":
    run_orchestrator()
