from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os


def get_amino_acid_count(file_path):
    """Counts the amino acids in an amino acid sequence from a fasta-file.

    Args:
        file_path: Path of the fasta-file to be read.

    Returns: A dictionary containing the amino acids as keys and the amino acid counts as values.

    """
    with open(file_path, 'r') as fasta_file:
        amino_acid_seq = ''
        for line in fasta_file:
            if not line.startswith('>'):
                amino_acid_seq += line.strip()
    return dict(Counter(amino_acid_seq))


def generate_output(count, file_name):
    """Saves the amino acid count as a csv-file and the graph of the amino acid distribution as a pdf-file.

    Args:
        count: A dictionary containing the amino acids as keys and the amino acid counts as values.
        file_name: Name of the fasta-file as string.

    Returns:

    """
    organism = file_name.replace('.fasta', '')
    count_df = pd.DataFrame({key: [val] for key, val in count.items()})
    count_df.to_csv(organism + '.csv')

    graph = plt.figure()
    plt.xlabel('Amino acides', size=11)
    plt.ylabel('Counts', size=11)
    plt.title(organism, size=12)
    sns.barplot(data=count_df)
    plt.show()
    graph.savefig(organism + '.pdf')


if __name__ == '__main__':
    for file in os.listdir('data'):
        aac = get_amino_acid_count('data/' + file)
        generate_output(aac, file)
