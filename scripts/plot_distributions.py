"""
Generates a kernel density estimate (KDE) plot for side-chain orientation 
angles within helical tripeptides, categorized by residue size.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Force non-interactive backend for headless environments
matplotlib.use("Agg")

def generate_density_visualization():
    # Handle command-line arguments or use defaults
    data_source = sys.argv[1] if len(sys.argv) > 1 else "final/angles.tsv"
    export_path = sys.argv[2] if len(sys.argv) > 2 else "final/ss_profile_HHH_for_arg_with_valid_runs.png"

    # Define strict ordering and color palette to match existing charts
    CATEGORIES = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]
    PALETTE = {
        "Tiny": "#e6e2d3",
        "Small": "#f2c879",
        "Intermediate": "#f28e5c",
        "Large": "#e4572e",
        "Bulky": "#b11226",
    }

    # Load and clean dataset
    dataset = pd.read_csv(data_source, sep="\t")

    # Normalize angles to the [-180, 180] range
    dataset["angle"] = ((dataset["angle"] + 180) % 360) - 180
    dataset = dataset[dataset["size_class"].isin(CATEGORIES)]

    total_samples = len(dataset)

    # Initialize the figure with specific styling for consistency
    fig, axis = plt.subplots(figsize=(10, 6))
    
    # Apply background colors verbatim from original
    fig.set_facecolor("#bdbdbd")
    axis.set_facecolor("#bdbdbd")

    # Define the X-axis resolution
    sampling_grid = np.linspace(-180, 180, 1000)

    # Calculate and plot KDE for each group
    for group in CATEGORIES:
        group_data = dataset[dataset["size_class"] == group]["angle"].values
        
        if len(group_data) < 2:
            continue
            
        # Use Silverman's rule for bandwidth selection to match previous output
        density_estimator = gaussian_kde(group_data, bw_method="silverman")
        y_values = density_estimator(sampling_grid)
        
        axis.plot(sampling_grid, y_values, color=PALETTE[group], linewidth=1.6, label=group)

    # Axis and Grid Formatting
    axis.set_xlim(-180, 180)
    grid_points = np.arange(-180, 181, 50)
    axis.set_xticks(grid_points)

    for point in grid_points:
        axis.axvline(point, color="blue", linestyle=":", linewidth=0.5, alpha=0.7)

    # Labels and Titles (using LaTeX formatting for special characters)
    axis.set_xlabel(r"Angle between adjacent C-$\alpha$ $\to$ Centroid vectors [°]")
    axis.set_ylabel("Norm. Freq. [A.U.]")
    axis.set_title(f"Tripeptide (XRX) in Helix (n = {total_samples})")

    # Legend Styling
    legend_box = axis.legend(loc="upper left", facecolor="white", edgecolor="black")
    legend_box.get_frame().set_linewidth(0.8)

    # Aesthetic cleanup
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    # File Export
    save_dir = os.path.dirname(export_path) or "."
    os.makedirs(save_dir, exist_ok=True)
    
    plt.tight_layout()
    plt.savefig(export_path, dpi=300, facecolor=fig.get_facecolor())
    
    print(f"Visualization saved to: {export_path}")
    print(f"Total count: {total_samples}")

if __name__ == "__main__":
    generate_density_visualization()
