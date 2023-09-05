import argparse
import pandas as pd

def replace_beta_with_entropy(pdb_file, csv_file):
    # Load CSV data
    csv_data = pd.read_csv(csv_file)
    csv_data['Residue number'] = csv_data['Residue number'].astype(int)
    entropy_dict = dict(zip(csv_data['Residue number'], csv_data['Entropy of residue']))
    print(entropy_dict)

    # Parse and modify PDB file
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()

    new_pdb_lines = []

    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            residue_number = int(line[22:26].strip())

            # Check if residue_number exists in CSV
            if residue_number in entropy_dict:
                entropy = entropy_dict[residue_number]
            else:
                entropy = 0.00  # default to zero if not found in csv

            # Replace beta factor and format it
            new_line = line[:60] + f"{entropy:>6.2f}" + line[66:]
            new_pdb_lines.append(new_line)
        else:
            new_pdb_lines.append(line)

    # Write the modified PDB lines to an output file
    with open('beta_' + pdb_file, 'w') as f:
        f.writelines(new_pdb_lines)

def main():
    parser = argparse.ArgumentParser(description="Replace beta column in PDB file with entropy values from a CSV file.")
    parser.add_argument("pdb_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("csv_file", type=str, help="Path to the input CSV file containing entropy values.")

    args = parser.parse_args()

    replace_beta_with_entropy(args.pdb_file, args.csv_file)

if __name__ == "__main__":
    main()
