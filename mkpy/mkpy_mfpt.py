import pandas as pd
import numpy as np
import argparse
import glob

def compute_distribution(data, monitored_cols):
    """
    Compute the mean and standard deviation of the distribution formed by the specified columns.
    """
    # Exclude all monitored columns from the distribution columns
    distribution_cols = [col for col in data.columns if col not in monitored_cols]
    
    distribution_data = data[distribution_cols].values.flatten()
    mean = np.nanmean(distribution_data)
    std = np.nanstd(distribution_data)
    
    return mean, std

def compute_first_passage_time(data, monitored_col, mean, std, sigma_range=1):
    """
    Compute the first passage time for the monitored column.
    """
    monitored_data = data[monitored_col]
    lower_bound = mean - sigma_range * std
    upper_bound = mean + sigma_range * std
    
    passage_time = next((index for index, value in monitored_data.iteritems() if lower_bound <= value <= upper_bound), None)
    
    return passage_time

def main(file_pattern, monitored_cols, sigma_range=1):
    """
    Main function to compute MFPT over multiple CSV files.
    """
    for monitored_col in monitored_cols:
        first_passage_times = []

        total_num_files = len(glob.glob(file_pattern))
        for csv_path in glob.glob(file_pattern):
            data = pd.read_csv(csv_path, header=None)
            
            mean, std = compute_distribution(data, monitored_cols)
            fpt = compute_first_passage_time(data, monitored_col, mean, std, sigma_range)
            
            if fpt is not None:
                first_passage_times.append(fpt)
        
        if not first_passage_times:
            print(f"No first passage times found for column {monitored_col}.")
            continue
        
        mfpt = np.nanmean(first_passage_times)
        num_files = len(first_passage_times)
        print(f"Column {monitored_col}: Averaged MFPT over {num_files}/{total_num_files} CSV files: {mfpt:.2f} frames")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute MFPT over multiple CSV files.')
    parser.add_argument('file_pattern', type=str, help='Pattern to match CSV files (e.g., "*/phi.csv").')
    parser.add_argument('--monitored_cols', nargs='+', type=int, default=[5], help='Indices of columns to monitor. Multiple columns can be specified.')
    parser.add_argument('--sigma_range', type=float, default=1, help='Range in terms of sigma to consider for the first passage.')
    args = parser.parse_args()

    main(args.file_pattern, args.monitored_cols, args.sigma_range)
