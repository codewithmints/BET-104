"""
Protein Structure Analysis Pipeline
Orchestrates the workflow from raw PDB files to finalized angle distribution plots.
"""

configfile: "config.yaml"

import glob
import os

# Configuration and Constants
TARGET_RESIDUE = config["target_aa"]
STRIDE_BINARY = config.get("stride_path", "/usr/local/bin/stride")

# Identify all available PDB datasets
PDB_IDENTIFIERS = [
    os.path.basename(f).replace(".pdb.gz", "")
    for f in glob.glob("pdbs/*.pdb.gz")
]

rule all:
    input:
        "results/ss_profile_HHH_for_arg_with_valid_runs.png",
        "results/pdbs.txt"

rule decompress_structure:
    input:
        compressed = "pdbs/{pdb_id}.pdb.gz"
    output:
        raw_pdb = temp("unzipped_pdbs/{pdb_id}.pdb")
    shell:
        "zcat {input.compressed} > {output.raw_pdb}"

rule assign_secondary_structure:
    input:
        raw_pdb = "unzipped_pdbs/{pdb_id}.pdb"
    output:
        stride_log = "stride_out/{pdb_id}.ss.out"
    params:
        bin = STRIDE_BINARY
    shell:
        "{params.bin} {input.raw_pdb} > {output.stride_log} 2>/dev/null || touch {output.stride_log}"

rule generate_triplet_contexts:
    input:
        stride_log = "stride_out/{pdb_id}.ss.out"
    output:
        context_data = "contexts/context_for_{aa}_in_{pdb_id}.tsv"
    params:
        residue = "{aa}"
    shell:
        "python scripts/extract_triplets.py {input.stride_log} {output.context_data} {params.residue}"

rule compute_geometric_angles:
    input:
        # Collects all context files before running the batch calculation
        dataset = expand("contexts/context_for_{aa}_in_{pdb_id}.tsv", pdb_id=PDB_IDENTIFIERS, aa=[TARGET_RESIDUE])
    output:
        tsv_out = "results/angles.tsv",
        log_out = "results/pdbs.txt"
    params:
        query_aa = TARGET_RESIDUE
    shell:
        "python scripts/angle_calc.py contexts {output.tsv_out} {params.query_aa} pdbs"

rule visualize_distribution:
    input:
        data_file = "results/angles.tsv"
    output:
        plot_image = "results/ss_profile_HHH_for_arg_with_valid_runs.png"
    shell:
        "python scripts/plot_distributions.py {input.data_file} {output.plot_image}"
