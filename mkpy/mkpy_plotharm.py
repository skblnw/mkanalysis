import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.optimize import curve_fit
import numpy as np
import os

# Define the sigmoid function
def sigmoid(x, L, k, x0):
    return L / (1 + np.exp(-k * (x - x0)))

def reversed_sigmoid(x, L, k, x0):
    return L - (L / (1 + np.exp(-k * (x - x0))))

def determine_direction_and_fit(data, column_to_fit, boundary=1000):
    """
    Determine the direction of the data and fit using appropriate sigmoid function.

    Parameters:
    - data: DataFrame containing the dataset.
    - column_to_fit: The column name to fit. Default is "success_rate".
    - boundary: Boundary for imaginary data points. Default is 1000.

    Returns:
    - Fitted parameters for the sigmoid function.
    - Modified data DataFrame.
    - Fitted data (smooth curve).
    """
    
    # Determine direction based solely on start and end data points
    direction = 'increasing' if data[column_to_fit].iloc[0] < data[column_to_fit].iloc[-1] else 'decreasing'

    # Choose the appropriate sigmoid function and initial parameters based on direction
    if direction == 'increasing':
        sigmoid_func = sigmoid
        initial_params = [1, 0.01, 100]

        # Add imaginary data points to the dataset
        data = data.append({"harm": 0, column_to_fit: 0}, ignore_index=True)
        data = data.append({"harm": boundary, column_to_fit: 1}, ignore_index=True)
        data.sort_values(by="harm", inplace=True)
    else:
        sigmoid_func = reversed_sigmoid
        initial_params = [1, -0.01, 100]

        # Add imaginary data points to the dataset
        data = data.append({"harm": 0, column_to_fit: 1}, ignore_index=True)
        data = data.append({"harm": boundary, column_to_fit: 0}, ignore_index=True)
        data.sort_values(by="harm", inplace=True)

    # Fit the data
    popt, pcov = curve_fit(sigmoid_func, data['harm'], data[column_to_fit], p0=initial_params)

    x_smooth = np.linspace(data['harm'].min(), data['harm'].max(), 1000)
    y_smooth = sigmoid_func(x_smooth, *popt)
    interpolation = {"harm": x_smooth, "sigmoid_fit": y_smooth}

    # Remove the imaginary data points
    data = data[data["harm"] != 0]
    data = data[data["harm"] != boundary]

    # Predict using the sigmoid function with the optimized parameters
    data['sigmoid_fit'] = sigmoid_func(data['harm'], *popt)

    return popt, pcov, data, interpolation

def fit_sigmoid(data, column_to_fit):
    # Fit the sigmoid function to the column_to_fit data
    popt, pcov, data, interpolation = determine_direction_and_fit(data, column_to_fit)
    errors_fit = np.sqrt(np.diag(pcov))

    # Extract the x0 value from the optimized parameters
    x0 = popt[2]
    error = errors_fit[2]
    print(f"Stiffness: {x0:.2f} ± {error:.2f}")
    # Determine the y-coordinate on the sigmoid curve at x0
    y_at_x0 = sigmoid(x0, *popt)
    # print(y_at_x0)

    return data, interpolation, x0, y_at_x0

def normalize_series(series):
    """
    Normalize a pandas Series to the range 0-1.
    """
    min_val = series.min()
    max_val = series.max()
    return (series - min_val) / (max_val - min_val)

def plot_primary_axis_with_fit(ax, data, show_interpolation):
    data, interpolation, x0, y_at_x0 = fit_sigmoid(data, 'success_rate')
    # ax.plot(data['harm'], data['sigmoid_fit'], 'b--', label='Sigmoid Fit')
    ax.axvline(x=x0, ymin=0, ymax=y_at_x0, color='k', linestyle=':', label=f'x0={x0:.2f}')
    if show_interpolation:
        ax.plot(interpolation['harm'], interpolation['sigmoid_fit'], 'k--', label='Sigmoid Fit Full')
        ax.set_xlim(0, data['harm'].max()*1.5)

def plot_secondary_axis_with_fit(ax, data, show_interpolation):
    data['convergence_norm'] = normalize_series(data['convergence'])
    ax.plot(data['harm'], data['convergence_norm'], color='tab:red', marker='o', linestyle='', label='Convergence')
    data, interpolation, x0, y_at_x0 = fit_sigmoid(data, 'convergence_norm')
    # ax.plot(data['harm'], data['sigmoid_fit'], 'r--', label='Sigmoid Fit')
    ax.axvline(x=x0, ymin=0, ymax=y_at_x0, color='k', linestyle=':', label=f'x0={x0:.2f}')
    if show_interpolation:
        ax.plot(interpolation['harm'], interpolation['sigmoid_fit'], 'k--', label='Sigmoid Fit Full')
        ax.set_xlim(0, data['harm'].max()*1.5)

def plot_data(args, data, save_path=None):
    fit_axes = args.fit_axes
    tick_font_size = 32

    fig, ax1 = plt.subplots(figsize=(7,8))
    
    if fit_axes in ['primary', 'both']:
        # ax1.set_ylabel('Convergence Time', color='tab:red')
        plot_secondary_axis_with_fit(ax1, data, args.show_interpolation)
        ax1.set_ylim(0, 1.2)
    else:
        ax1.set_ylabel('Convergence Time (ns)', color='tab:red')
        ax1.errorbar(data['harm'], data['convergence']/10, yerr=data['error'], color='tab:red', marker='o', linestyle='-', label='Convergence')
        ax1.set_ylim(0, 100)
    ax1.set_xlim(0, 190)
    ax1.set_xticks([0,50,100,150])
    ax1.tick_params(axis='y', labelcolor='tab:red', labelsize=tick_font_size)
    ax1.tick_params(axis='x', labelsize=tick_font_size)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    
    ax2 = ax1.twinx()
    ax2.set_zorder(1)  # Set zorder of ax2 to be behind ax1
    ax1.set_zorder(2)  # Set zorder of ax1 to be in front of ax2
    ax1.patch.set_visible(False)  # Make ax1 background transparent

    # ax2.set_ylabel('Success Rate', color='tab:blue')
    ax2.fill_between(data['harm'], 0, data['success_rate'], color='tab:blue', alpha=0.3)
    # ax2.bar(data['harm'], data['success_rate'], color='tab:blue', width=4.5, alpha=0.6, label='Success Rate')
    if fit_axes in ['secondary', 'both']:
        plot_primary_axis_with_fit(ax2, data, args.show_interpolation)
    ax2.set_ylim(0, 1.2)
    ax2.tick_params(axis='y', labelcolor='tab:blue', labelsize=tick_font_size)
    ax2.spines['top'].set_visible(False)
    
    plt.tight_layout()
    # plt.show()
    plt.savefig(save_path, format='pdf', bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser(description="Plot data from CSV file.")
    parser.add_argument("csv_file", help="Path to the CSV file.")
    parser.add_argument("--fit-axes", default="both", help="Which axis to perform curve fitting.")
    parser.add_argument("--show-interpolation", action='store_true', help='If set, an interpolated line will be shown.')
    parser.add_argument("--pdf-file", help="Path to save the PDF file. By default, it uses the prefix of the input filename.")
    args = parser.parse_args()

    # Load the data from the CSV file
    data = pd.read_csv(args.csv_file)
    # Calculate success rate
    data['success_rate'] = data['successful'] / data['trials']

    # Determine the PDF filename
    if args.pdf_file:
        pdf_filename = args.pdf_file
    else:
        # Default PDF filename using the prefix of the input CSV filename
        pdf_filename = os.path.splitext(args.csv_file)[0] + ".pdf"

    # Plot the data and save as PDF
    plot_data(args, data, save_path=pdf_filename)

if __name__ == "__main__":
    main()
