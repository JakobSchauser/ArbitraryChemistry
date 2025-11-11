import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
import csv

def main(L, D):
    # Get script directory
    script_dir = Path(__file__).parent
    
    # Load data from outputs in script directory
    data_dir = script_dir / "outputs" / "latticeTimeseries" / f"L_{L}_D_{D}"
    
    # Find all lattice files
    lattice_files = list(data_dir.glob("lattice_*.txt"))
    
    # Collect data: list of (timestep, global_counts)
    data_list = []
    for file_path in lattice_files:
        # Extract timestep from filename (e.g., lattice_1000.txt -> 1000)
        t = int(file_path.stem.split('_')[1])
        
        # Read and parse lattice file (new CSV format)
        with open(file_path, 'r') as f:
            # Skip comment lines starting with '#'
            lines = [line for line in f if not line.startswith('#')]
        
        # Initialize sum_counts for chemicals 0 to 100
        sum_counts = np.zeros(101, dtype=int)  # 0 to 100
        
        # Parse CSV data (skip header line)
        reader = csv.DictReader(lines)
        for row in reader:
            chemical = int(row['chemical'])
            count = int(row['count'])
            sum_counts[chemical] += count
        
        data_list.append((t, sum_counts))
    
    # Sort by timestep
    data_list.sort(key=lambda x: x[0])
    
    # Extract steps and counts
    steps = np.array([item[0] for item in data_list])
    counts = np.array([item[1] for item in data_list])
    
    # counts[:, 0] is empty, counts[:, 1:] are chemicals 1 to 100
    n_chemicals = counts.shape[1] - 1  # Exclude empty (chemical 0)
    
    # Find the highest chemical index with non-zero counts
    max_chemical_index = 0
    for c in range(1, n_chemicals + 1):
        if np.max(counts[:, c]) > 0:
            max_chemical_index = c
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Rainbow colormap: red for heavy, blue/violet for light, up to max_chemical_index
    colors = cm.rainbow(np.linspace(0, 1, max_chemical_index + 1))
    
    # Plot each chemical up to max_chemical_index
    for c in range(1, max_chemical_index + 1):
        ax.plot(steps, counts[:, c], label=f'Chemical {c}', color=colors[c], linewidth=1.5)
    
    # Plot empty separately with black
    ax.plot(steps, counts[:, 0], label='Empty (0)', color='black', linewidth=1.5, linestyle='--')
    
    ax.set_xlabel('Timestep', fontsize=14)
    ax.set_ylabel('Molecule Count', fontsize=14)
    ax.set_yscale('log')
    ax.set_title(f'Chemical Populations Over Time (L={L}, D={D})', fontsize=16)
    ax.grid(True, alpha=0.3)
    
    # Save figure to plots in script directory
    output_dir = script_dir / "plots" / "latticeTimeseries"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"populations_L_{L}_D_{D}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")
    plt.close()

if __name__ == "__main__":
    L = 25  # Lattice size
    D = 1   # Diffusion coefficient

    main(L, D)