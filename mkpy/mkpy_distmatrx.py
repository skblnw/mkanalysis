"""
# PDB Distance Matrix Calculator

This script takes an input PDB file and calculates two distance matrices for CA and CB atoms based on specified residue ranges for each chain.

## Features

1. Reads PDB file and extracts coordinates for CA and CB atoms within the specified residue ranges for each chain.
2. Computes the distance matrix for CA and CB atoms separately.
3. Writes the distance matrices to files.
4. Plots the distance matrices as images with a grayscale color map.

## Usage

To run the script, use the following command:

```
python script.py <input_PDB_file>
```

The script will create four output files:

- ca_distance_matrix.txt: Distance matrix for CA atoms.
- cb_distance_matrix.txt: Distance matrix for CB atoms.
- ca_distance_matrix.png: Plot of the distance matrix for CA atoms.
- cb_distance_matrix.png: Plot of the distance matrix for CB atoms.

## Customization

To change the residue ranges for each chain, modify the PARAMETERS dictionary in the script. The key is the chain identifier (e.g., "A", "C"), and the value is a list of two integers representing the start and end residue numbers (inclusive).

Example:

```
PARAMETERS = {
    "A": [1, 180],
    "C": [1, 9]
}
```

## Requirements

This script is written in Python and requires Python 3 and the following libraries to run:

- math
- matplotlib

"""

import sys
import math
import matplotlib.pyplot as plt

# Define the parameters in the script
PARAMETERS = {
    "A": [1, 180],
    "C": [1, 9]
}

def read_pdb(filename, parameters):
    with open(filename, 'r') as file:
        lines = file.readlines()
    ca_atoms = []
    cb_atoms = []
    for line in lines:
        if line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            chain_id = line[21].strip()
            res_num = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            if chain_id in parameters:
                if parameters[chain_id][0] <= res_num <= parameters[chain_id][1]:
                    if atom_name == "CA":
                        ca_atoms.append((x, y, z))
                    elif atom_name == "CB":
                        cb_atoms.append((x, y, z))
    return ca_atoms, cb_atoms

def euclidean_distance(coord1, coord2):
    return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)

def compute_distance_matrix(atoms):
    n_atoms = len(atoms)
    matrix = [[0.0 for _ in range(n_atoms)] for _ in range(n_atoms)]

    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            distance = euclidean_distance(atoms[i], atoms[j])
            matrix[i][j] = distance
            matrix[j][i] = distance

    return matrix

def write_distance_matrix(filename, matrix):
    with open(filename, 'w') as file:
        for row in matrix:
            file.write(" ".join("{:.2f}".format(x) for x in row) + "\n")

def plot_distance_matrix(matrix, title, output_file):
    if not matrix or not matrix[0]:
        print(f"Warning: Empty matrix for {title}. Skipping plot generation.")
        return

    plt.imshow(matrix, cmap='Greys', interpolation='nearest')
    plt.colorbar()
    plt.title(title)
    plt.savefig(output_file)
    plt.clf()

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <PDB filename>")
        exit(1)

    pdb_filename = sys.argv[1]
    ca_atoms, cb_atoms = read_pdb(pdb_filename, PARAMETERS)

    ca_distance_matrix = compute_distance_matrix(ca_atoms)
    cb_distance_matrix = compute_distance_matrix(cb_atoms)

    write_distance_matrix("ca_distance_matrix.txt", ca_distance_matrix)
    write_distance_matrix("cb_distance_matrix.txt", cb_distance_matrix)

    plot_distance_matrix(ca_distance_matrix, "CA Distance Matrix", "ca_distance_matrix.png")
    plot_distance_matrix(cb_distance_matrix, "CB Distance Matrix", "cb_distance_matrix.png")

if __name__ == "__main__":
    main()
