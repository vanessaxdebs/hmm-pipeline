import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

def calculate_entropy(alignment_file):
    with open(alignment_file, 'r') as f:
        sequences = [line.strip() for line in f.readlines() if not line.startswith('>')]

    seq_length = len(sequences[0])
    entropy_values = []

    for i in range(seq_length):
        column = [seq[i] for seq in sequences]
        count = Counter(column)
        total = len(column)
        prob = [count[char] / total for char in count]
        entropy = -sum(p * np.log2(p) for p in prob)
        entropy_values.append(entropy)

    return entropy_values

def plot_entropy(entropy_values, output_file):
    plt.plot(entropy_values)
    plt.xlabel('Position in Alignment')
    plt.ylabel('Entropy')
    plt.title('Entropy Plot (Per-Column)')
    plt.savefig(output_file)
    plt.show()

# Usage
entropy_values = calculate_entropy("example/msa.aln")
plot_entropy(entropy_values, "entropy_plot.png")

