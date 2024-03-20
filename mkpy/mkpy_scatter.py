import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np

def compute_alpha(row_idx, max_row_idx, alpha_min=0.1, alpha_max=1.0):
    """Computes the alpha (transparency) value based on row index."""
    normalized = row_idx / max_row_idx
    return alpha_min + (alpha_max - alpha_min) * normalized

def plot_scatter(file1, file2, indices=None, colormap='viridis'):
    # Load data
    with open(file1, 'r') as f:
        phi_angles = [list(map(float, line.strip().split(','))) for line in f]
        
    with open(file2, 'r') as f:
        psi_angles = [list(map(float, line.strip().split(','))) for line in f]
    
    # If indices is None, select all columns.
    if indices is None:
        indices = list(range(len(phi_angles[0])))
    
    # Generate color list
    cmap = plt.get_cmap(colormap)
    colors = [cmap(i) for i in np.linspace(0, 1, len(indices))]

    plt.figure(figsize=(10, 8))

    # Select specified indices and scatter plot with corresponding color
    for ii, (idx, color) in enumerate(zip(indices, colors)):
        phi_column = [360-angles[idx] for angles in phi_angles]
        psi_column = [360-angles[idx] for angles in psi_angles]
        plt.scatter(phi_column, psi_column, color=color, label=f"Pos {ii+1}", alpha=0.3, edgecolors='none', s=60)
        # max_row_idx = len(phi_angles) - 1  # Assuming phi_angles and psi_angles have the same length
        # alphas = [compute_alpha(row_idx, max_row_idx) for row_idx in range(len(phi_column))]
        # for phi, psi, alpha in zip(phi_column, psi_column, alphas):
        #     plt.scatter(phi, psi, color=color, alpha=alpha, label=f"Index {idx}" if phi == phi_column[0] else "")

    plt.legend(fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=24)  # For major ticks
    # plt.tick_params(axis='both', which='minor', labelsize=12)  # For minor ticks (if any)
    # plt.xticks(ticks=np.arange(0,420,60), labels=['0°', '60°', '120°', '180°', '240°', '300°', '360°'])
    # plt.yticks(ticks=np.arange(0,420,60), labels=['0°', '60°', '120°', '180°', '240°', '300°', '360°'])
    plt.xticks(ticks=np.arange(0,420,60), labels=['-180°', '-120°', '-60°', '0°', '60°', '120°', '180°'])
    plt.yticks(ticks=np.arange(0,420,60), labels=['-180°', '-120°', '-60°', '0°', '60°', '120°', '180°'])
    # plt.title('Scatter plot of Phi vs Psi angles')
    # plt.xlabel('Phi Angles')
    # plt.ylabel('Psi Angles')
    # plt.grid(True)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scatter plot of Phi and Psi angles over a molecular trajectory.")
    parser.add_argument("file1", type=str, help="File containing Phi angles.")
    parser.add_argument("file2", type=str, help="File containing Psi angles.")
    parser.add_argument("--index", type=str, help="Indices to plot, separated by ',' or a range using '-'. E.g., 5,6,7 or 5-10.", default=None)
    parser.add_argument("--lastN", type=int, help="Plot the last N columns.", default=None)
    parser.add_argument("--colormap", type=str, default="tab10", help="Colormap to use for differentiating columns. Default is 'viridis'.")


    args = parser.parse_args()
    
    # Determine columns to be plotted based on provided arguments
    if args.index:
        if '-' in args.index:
            start, end = map(int, args.index.split('-'))
            indices = list(range(start, end+1))
        else:
            indices = list(map(int, args.index.split(',')))
    elif args.lastN:
        # Here we will assume both files have the same number of columns for simplicity.
        # If not, additional logic might be required.
        with open(args.file1, 'r') as f:
            columns = len(f.readline().strip().split(','))
        indices = list(range(columns - args.lastN, columns))
    else:
        indices = None
    
    print(indices)
    plot_scatter(args.file1, args.file2, indices, args.colormap)
