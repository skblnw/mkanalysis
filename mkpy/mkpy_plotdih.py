import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import argparse

def moving_average(data, window_size=5):
    """
    Calculate the moving average of the given list or series.
    """
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def plot_dih_angles(csv_file, save_as_tiff=False, tiff_filename="output.tiff", show_fitting=True):
    """
    Plot dih angles from a CSV file.
    """
    data = pd.read_csv(csv_file, header=None)

    fig, ax = plt.subplots(1, 2, figsize=(15, 6))

    cmap = plt.get_cmap("tab20")

    for idx, col in enumerate(data.columns):
        ax[0].plot(data[col], alpha=0.1, color=cmap(idx))
        ax[0].plot(moving_average(data[col]), label=f"dih {idx}", color=cmap(idx))

    ax[0].set_xlim([0, len(data) * 1.2])  # Extend x-axis by 20%
    ax[0].legend(loc='upper right')
    ax[0].set_title('Profiles with Moving Average')
    ax[0].set_xlabel('Number of Frames')
    ax[0].set_ylabel('Value')
    ax[0].set_ylim([0, 360])

    for idx, col in enumerate(data.columns):
        column_data = data[col].dropna()
        
        if show_fitting:
            mu, std = norm.fit(column_data)
            xmin, xmax = plt.xlim()
            y = np.linspace(0, 360, 100)
            p = norm.pdf(y, mu, std)
            ax[1].plot(p * (xmax-xmin) + xmin, y, 'k', linewidth=2)

        ax[1].hist(column_data, bins=50, density=True, alpha=0.6, range=(0, 360), label=f"dih {idx}", orientation='horizontal', color=cmap(idx))

    ax[1].set_title('Histograms with Fitted Distribution')
    ax[1].set_ylabel('Value')
    ax[1].set_xlabel('Density')
    ax[1].set_ylim([0, 360])
    fig.suptitle(f"{data.shape[1]} dih angles", fontsize=16)

    plt.tight_layout()

    if save_as_tiff:
        plt.savefig(tiff_filename, format='tiff', dpi=300)

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot dih angles from a CSV file.')
    parser.add_argument('csv_file', type=str, help='Path to the input CSV file.')
    parser.add_argument('--output', type=str, default=None, help='Name of the TIFF file to save the plot.')
    parser.add_argument('--no-fitting', action='store_true', help='Do not show normal fitting on histograms.')
    args = parser.parse_args()

    plot_dih_angles(args.csv_file, save_as_tiff=bool(args.output), tiff_filename=args.output, show_fitting=not args.no_fitting)
