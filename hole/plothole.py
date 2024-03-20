import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_data(args):
    # Define colors and font sizes
    # line_color = '#1f77b4'  # Example blue color
    # line_color = '#63CF8A'  # Example blue color
    # line_color = '#AD71C1'  # Example blue color
    # line_color = '#DFDE53'  # Example blue color
    line_color = '#F2658B'  # Example blue color
    axis_font_size = 12
    legend_font_size = 10
    title_font_size = 14
    x_ticks_font_size = 18  # Font size for x-axis ticks

    # Load the data
    data = pd.read_csv(args.file_path, header=None, names=['y', 'x'], delimiter=' ')

    # Create the third column by negating the second column
    data['x_mirrored'] = -data['x']

    # Plotting
    plt.figure(figsize=(2.5, 8))

    # Plot original and mirrored data
    plt.plot(data['x'], data['y'], label='Original Line', color=line_color)
    plt.plot(data['x_mirrored'], data['y'], label='Mirrored Line', color=line_color)

    # Add shaded area between the two lines
    plt.fill_betweenx(data['y'], data['x'], data['x_mirrored'], color=line_color, alpha=0.2)

    # Axis labels, title, and ticks formatting
    # plt.xlabel('X axis', fontsize=axis_font_size)
    # plt.ylabel('Y axis', fontsize=axis_font_size)
    # plt.title('Plot of Y vs. X and its Mirror', fontsize=title_font_size)
    plt.gca().yaxis.set_visible(False)  # Hide y-axis ticks
    plt.gca().tick_params(axis='x', labelsize=x_ticks_font_size)  # Font size for x-axis ticks

    # Hide top, left, and right spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # Legend
    # plt.legend(fontsize=legend_font_size)

    # Optional: Vertical line at x=0 for reference
    plt.axvline(x=0, color='grey', lw=1, linestyle='--')

    # Decision on showing or saving the figure
    if args.output:
    	outfile=args.output
    	plt.savefig(outfile, bbox_inches='tight', pad_inches=.02)

    if args.interactive:
    	plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot data from a file with mirrored x-axis values.')
    parser.add_argument('file_path', type=str, help='Path to the file containing the data to be plotted')

    io_group = parser.add_mutually_exclusive_group(required=True)
    io_group.add_argument('-o', '--output', type=str, help='PDF output file')
    io_group.add_argument('-i', '--interactive', action='store_true', help='Launches an interactive matplotlib session')

    args = parser.parse_args()
    
    plot_data(args)
