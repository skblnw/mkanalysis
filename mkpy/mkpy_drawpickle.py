import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def visualize_matrices_from_pickle(pickle_file_path):
    # Load dictionary from pickle file
    with open(pickle_file_path, 'rb') as f:
        matrix_dict = pickle.load(f)
    
    for keyword, matrix in matrix_dict.items():
        # Plotting the heatmap for each matrix
        plt.figure(figsize=(10,8))
        sns.heatmap(matrix, cmap='viridis', annot=False)  # 'annot=True' will annotate each cell with its value
        plt.title(f'Heatmap for {keyword}')
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Visualize matrices from a pickle file using a heatmap.")
    parser.add_argument("pickle_file_path", type=str, help="Path to the pickle file containing matrices.")

    args = parser.parse_args()
    
    visualize_matrices_from_pickle(args.pickle_file_path)

if __name__ == "__main__":
    main()
