import pandas as pd
from itertools import combinations
import argparse
import glob

def compute_correlations_and_merge(filepaths):
    """
    Compute correlation coefficients between multiple time series data files and merge them into a new csv.
    
    Parameters:
    - filepaths : list of str
        Paths to the time series data files.
        
    Returns:
    - dict
        Dictionary of correlation coefficients with file pairs as keys.
    """
    
    # Custom parser to handle lines starting with '#'
    def custom_parser(filepath, col_name):
        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f.readlines() if not line.startswith("#")]
        data = [tuple(map(float, line.split())) for line in lines]
        return pd.DataFrame(data, columns=["time", col_name])
    
    # 1. Read the time series data files using the custom parser
    dfs = [custom_parser(filepath, filepath) for filepath in filepaths]
    
    # 2. Merge them based on a common index (assuming 'time')
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on="time", how="inner")
    
    # 3. Compute the correlation coefficient for each pair of files
    correlation_results = {}
    for file1, file2 in combinations(filepaths, 2):
        correlation = merged_df[file1].corr(merged_df[file2])
        correlation_results[(file1, file2)] = correlation
    
    # 4. Write the merged data into a new CSV with filenames as the column names
    # merged_df.to_csv("merged_data.csv", index=False)
    
    return correlation_results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute correlation coefficients between multiple time series data files and merge them.")
    parser.add_argument("pattern", type=str, help="Wildcard pattern to specify input files (e.g., 'data*.txt')")

    args = parser.parse_args()
    
    # Use glob to get file paths that match the wildcard pattern
    filepaths = glob.glob(args.pattern)

    filepaths = sorted(filepaths)

    res = []
    correlations = compute_correlations_and_merge(filepaths)
    for (file1, file2), corr in correlations.items():
        print(f"Correlation coefficient between {file1} and {file2}: {corr:.2}")
        res.append(round(corr,2))
    print(res)
