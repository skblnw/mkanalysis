import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def moving_average(data, window_size):
    """
    Compute moving average using a specified window size.
    """
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def plot_csv_data(csv_path, window_size=5, step_size=10):
    """
    Plot the data from the CSV file with moving average, boxplot of moving values, and cumulative boxplot.
    """
    data = pd.read_csv(csv_path, header=None)

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
    
    # Line Plot with Moving Average for each column
    ax = axes[0]
    for col in range(data.shape[1]):
        ax.plot(data.index, data[col], alpha=0.4, label=f'Dihedral {col}')
        avg_data = moving_average(data[col], window_size)
        ax.plot(range(window_size - 1, len(data)), avg_data, label=f'Moving Avg Dihedral {col}')
    ax.set_title('Dihedral Angle vs Frame')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Angle (degrees)')
    # ax.legend()

    # Boxplot of Moving Values (considering all columns)
    ax = axes[1]
    box_data = [data.iloc[i-window_size+1:i+window_size].values.flatten() for i in range(window_size, len(data), step_size)]
    ax.boxplot(box_data, vert=True)
    ax.set_title('Boxplot of Moving Values')
    ax.set_xlabel('Data Points (Every 10 Frames)')
    ax.set_ylabel('Angle (degrees)')
    
    # Cumulative Boxplot (considering all columns)
    ax = axes[2]
    cum_data = [data.iloc[i:].values.flatten() for i in range(window_size, len(data), step_size)]
    ax.boxplot(cum_data, vert=True)
    ax.set_title('Cumulative Boxplot')
    ax.set_xlabel('Data Points (Every 10 Frames)')
    ax.set_ylabel('Angle (degrees)')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot data from a CSV file with various visualizations.')
    parser.add_argument('csv_path', type=str, help='Path to the CSV file.')
    parser.add_argument('--window_size', type=int, default=5, help='Window size for the moving average.')
    parser.add_argument('--step_size', type=int, default=10, help='Step size for the boxplots.')
    args = parser.parse_args()

    plot_csv_data(args.csv_path, args.window_size, args.step_size)
