import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path

def main(N, S):
    # Get script directory
    script_dir = Path(__file__).parent
    
    # Load data from outputs in script directory
    data_dir = script_dir / "outputs" / "timeseries" / f"N_{N}_S_{S}"
    counts_file = data_dir / "counts.tsv"
    
    # Read counts file
    with open(counts_file, 'r') as f:
        lines = f.readlines()
    
    # Parse header and volatility
    header = lines[0].strip().split('\t')
    volatility_line = lines[1].strip().split('\t')
    
    # Parse data (skip header and volatility rows)
    data = []
    for line in lines[2:]:
        values = [float(x) for x in line.strip().split('\t')]
        data.append(values)
    
    data = np.array(data)
    steps = data[:, 0]
    counts = data[:, 1:]  # All chemical counts (including empty at index 0)
    
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
    # ax.set_yscale('log')
    ax.set_title(f'Chemical Populations Over Time (N={N}, S={S})', fontsize=16)
    ax.grid(True, alpha=0.3)
    
    # Save figure to plots in script directory
    output_dir = script_dir / "plots" / "timeseries"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"populations_N_{N}_S_{S}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")
    plt.close()

if __name__ == "__main__":
    N = 100000  # Number of molecules
    S = 1000   # Supply rate

    main(N, S)