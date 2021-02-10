from collections import Counter
import sys


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


if __name__ == '__main__':
    file = sys.argv[1][1:-1]
    # file = 'exercise_2/data/Homo_sapiens.fasta'
    aac = get_amino_acid_count(file)
    print(aac)
