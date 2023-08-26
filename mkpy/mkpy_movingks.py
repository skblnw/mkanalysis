import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import argparse
import matplotlib.pyplot as plt
import glob

# =============================
# Data Processing Functions
# =============================
def cumulative_ks_test(data, late_fraction=0.1):
    """
    Perform a cumulative Kolmogorov-Smirnov test.
    
    Parameters:
    - data: DataFrame containing the dataset.
    - late_fraction: Fraction of data considered as late-stage.
    
    Returns:
    - List of tuples with frame numbers and their corresponding p-values.
    """

    late_start = int(len(data) * (1 - late_fraction))
    # late_data = data.iloc[late_start:].values.flatten()
    late_data = bootstrap_optimal_distribution(data)
    results = []

    for i in range(late_start):
        cumulative_data = data.iloc[i:].values.flatten()
        _, p_value = ks_2samp(cumulative_data, late_data)
        results.append((i, p_value))
    
    return results

def bootstrap_optimal_distribution(data, late_fraction=0.1, n_iterations=1000):
    """
    Use bootstrap resampling to determine the optimal late-stage distribution with the smallest variance.

    Parameters:
    - data: DataFrame containing the data
    - late_fraction: Fraction of frames to be considered as late-stage
    - n_iterations: Number of bootstrap iterations

    Returns:
    - The optimal distribution with the smallest variance
    """
    late_start = int(len(data) * (1 - late_fraction))
    late_data = data.iloc[late_start:]
    
    optimal_distribution = None
    min_variance = float('inf')

    for _ in range(n_iterations):
        # Randomly sample columns
        sampled_data = late_data.sample(n=late_data.shape[1] - 4, replace=False, axis=1)
        combined_distribution = sampled_data.values.flatten()

        # Compute variance for this iteration
        current_variance = np.var(combined_distribution)

        # Update optimal distribution if current variance is smaller
        if current_variance < min_variance:
            min_variance = current_variance
            optimal_distribution = combined_distribution

    return optimal_distribution

def compute_outlier_percentage(data, late_fraction=0.1):
    """
    Compute the outlier percentage for each data point in a cumulative box plot.
    
    Parameters:
    - data: DataFrame containing the dataset.
    - late_fraction: Fraction of data considered as late-stage.
    
    Returns:
    - List of tuples with frame numbers and their corresponding outlier percentages.
    """

    # List to store results
    results = []

    # Iterate over data points
    for i in range(10, int(len(data) * (1 - late_fraction))):  # Leaving last data point to have meaningful analysis
        cumulative_data = data.iloc[i:].values.flatten()

        # Compute quartiles and interquartile range
        q1 = np.percentile(cumulative_data, 25)
        q3 = np.percentile(cumulative_data, 75)
        iqr = q3 - q1

        # Determine bounds for outliers
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr

        # Count outliers
        # outlier_count = np.sum((cumulative_data < lower_bound) | (cumulative_data > upper_bound))
        outlier_count = np.sum((cumulative_data < lower_bound))
        
        # Compute percentage of outliers
        percentage = (outlier_count / len(cumulative_data)) * 100
        results.append((i, percentage))
    
    return results

def exponential_smoothing(data, alpha=0.1, threshold=0.01):
    """
    Apply exponential smoothing to a dataset.
    
    Parameters:
    - data: List of data values.
    - alpha: Smoothing factor (between 0 and 1).
    
    Returns:
    - List of smoothed data values.
    """

    # Initialize smoothed values list
    smoothed = [data[0]]

    # Compute smoothed values
    for i in range(1, len(data)):
        smoothed_value = alpha * data[i] + (1 - alpha) * smoothed[i - 1]
        smoothed.append(smoothed_value)

    return smoothed

def determine_convergence_ks(frames, data, threshold_value, threshold_type='value'):
    """
    Determine the frame where data reaches a given threshold.
    
    Parameters:
    - frames: List of frame numbers.
    - data: List of data values (smoothed percentages or p-values).
    - threshold_value: Threshold value for determining convergence.
    - threshold_type: 'value' for direct threshold or 'percentage' for a percentage of total decrease.
    
    Returns:
    - Frame number at which the data reaches the specified threshold.
    """

    if threshold_type == 'value':
        for frame, value in zip(frames, data):
            if value > threshold_value:
                return frame
    return None

def determine_convergence_box(frames, data, threshold_value, threshold_type='value', activation_threshold=5, late_fraction=0.05):
    """
    Determine the frame where data reaches a given threshold.
    
    Parameters:
    - frames: List of frame numbers.
    - data: List of data values (smoothed percentages or p-values).
    - threshold_value: Threshold value for determining convergence.
    - threshold_type: 'value' for direct threshold or 'percentage' for a percentage of total decrease.
    
    Returns:
    - Frame number at which the data reaches the specified threshold.
    """

    if threshold_type == 'value':
        for frame, value in zip(frames, data):
            if value <= threshold_value:
                return frame
    elif threshold_type == 'percentage':

        # Check if the maximum values are below the activation threshold
        if np.max(data) < activation_threshold:
            return None

        # Check if the minimum values are above the activation threshold
        if np.min(data) > activation_threshold:
            return None

        # Check if the late-stage values are below the activation threshold
        # late_start = int(len(data) * (1 - late_fraction))
        # late_data = data[late_start:]
        # if np.mean(late_data) > activation_threshold:
        #     return None

        # decrease = data[0] - data[-1]
        decrease = np.max(data) - np.min(data)
        target_value = data[0] - (decrease * threshold_value / 100)
        for frame, value in zip(frames, data):
            if value <= target_value:
                return frame
    return None

def bootstrap_mean_uncertainty(data, n_iterations=100):
    """
    Compute the mean and uncertainty of the data using bootstrap resampling.
    
    Parameters:
    - data: List of non-None values.
    - n_iterations: Number of bootstrap iterations.
    
    Returns:
    - Tuple containing the mean and its uncertainty (standard deviation).
    """
    bootstrap_means = []
    for _ in range(n_iterations):
        # Resample the data with replacement
        resampled_data = np.random.choice(data, size=len(data), replace=True)
        bootstrap_means.append(np.mean(resampled_data))
    
    # Compute mean and standard deviation of the bootstrapped means
    mean_value = np.mean(bootstrap_means)
    uncertainty = np.std(bootstrap_means)
    
    return mean_value, uncertainty

# =============================
# Plotting Functions
# =============================

def plot_results(frames, data, smoothed_data=None, fitted_data=None, test_name=""):
    """
    Plot the results of the specified test.
    
    Parameters:
    - frames: List of frame numbers.
    - data: List of data values.
    - smoothed_data: List of smoothed data values.
    - test_name: Name of the test being plotted.
    """

    plt.figure(figsize=(10, 6))
    
    # Plot original data
    plt.plot(frames, data, marker='o', linestyle='-', label="Original Data")
    
    # Plot smoothed data, if provided
    if smoothed_data is not None:
        plt.plot(frames, smoothed_data, marker='o', linestyle='-', color='red', label="Smoothed Data")
    
    # Plot fitted data, if provided
    if fitted_data is not None:
        plt.plot(frames, fitted_data, linestyle='--', color='blue', label="Fitted Data")
    
    plt.xlabel('Frame')
    plt.ylabel('Value')
    plt.title(f'{test_name} Data')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# =============================
# Main Execution Logic
# =============================

# Error handling function
def validate_input_data(args):
    # Read the data
    data = pd.read_csv(args.csv_path, header=None)

    # Check if data is empty or not in expected format
    if data.empty:
        raise ValueError("The provided CSV file is empty.")

    # Check for valid late_fraction value
    if not 0 < args.late_fraction < 1:
        raise ValueError("The late_fraction should be between 0 and 1.")

    return data

def main(args):
    file_paths = glob.glob(args.csv_path_pattern)
    total_files = len(file_paths)
    successful_activations = 0
    convergence_frames = []

    for path in file_paths:
        args.csv_path = path
        try:
            data = validate_input_data(args)
        except Exception as e:
            print(f"Error processing file {path}: {e}")
            continue

        # Determine the test to perform
        if args.test == 'ks':
            results = cumulative_ks_test(data, args.late_fraction)
            test_name = "Cumulative KS Test"
            frames, p_values = zip(*results)
            smoothed = exponential_smoothing(p_values)
            convergence_frame = determine_convergence_ks(frames, p_values, args.threshold_value)
            # print(convergence_frame)
        elif args.test == 'box':
            results = compute_outlier_percentage(data, args.late_fraction)
            test_name = "Cumulative Box Test"
            frames, percentages = zip(*results)
            smoothed = exponential_smoothing(percentages)
            convergence_frame = determine_convergence_box(frames, smoothed, args.threshold_value, args.box_check_method, args.activation_threshold)
            # print(f"path: {path}, cf: {convergence_frame}")
        
        if not args.no_plotting:
            plot_results(frames, percentages if args.test == 'box' else p_values, smoothed_data=smoothed, test_name=test_name)

        if convergence_frame is not None:
            successful_activations += 1
            convergence_frames.append(convergence_frame)

    if convergence_frames:
        average_convergence, uncertainty = bootstrap_mean_uncertainty(convergence_frames)
        print(f"Averaged Convergence for {test_name} over {successful_activations}/{total_files} trials: {average_convergence:.2f} ± {uncertainty:.2f}")
    else:
        print("No convergence was found")

# Argument parsing is moved to __main__
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze data convergence using KS and Box tests.')
    parser.add_argument('csv_path_pattern', type=str, help='Pattern for the path of the CSV files to be processed.')
    parser.add_argument('--no-plotting', action='store_true', help='If set, no plots will be generated.')
    parser.add_argument('--test', type=str, choices=['ks', 'box'], default='ks', help='Test to perform: "ks" for Kolmogorov-Smirnov test and "box" for Box test. Default is "ks".')
    parser.add_argument('--late_fraction', type=float, default=0.1, help='Fraction of data considered as late-stage. Default is 0.1.')
    parser.add_argument('--box_check_method', type=str, choices=['abs', 'percentage'], default='percentage', help='Method to use for Box test activation checking: "abs" for absolute value threshold and "percentage" for percentage threshold. Default is "percentage".')
    parser.add_argument('--threshold_value', type=float, default=0.05, help='Threshold value for determining convergence in KS test. Default is 0.05.')
    parser.add_argument('--activation_threshold', type=float, default=10, help='Threshold below which the late-stage values must fall to consider the system activated for Box test when using percentage threshold type. Default is 5%.')
    
    args = parser.parse_args()
    main(args)
