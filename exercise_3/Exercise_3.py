import pandas as pd
import plotly.graph_objects as go
from collections import deque


def get_amino_acid_sequence(file_path):
    """Determines the amino acid sequence of a fasta-file.

    Args:
        file_path: Path to fasta-file as string.

    Returns: amino acid sequence in 1-letter code as string

    """
    with open(file_path, 'r') as fasta_file:
        amino_acid_seq = ''
        for line in fasta_file:
            if not line.startswith('>'):
                amino_acid_seq += line.strip()
    return amino_acid_seq


def get_amino_acid_properties(file_path):
    """Reads the hydropathy of amino acids from a csv-file.

    Args:
        file_path: Path to the csv-file as string.

    Returns: hydropathy of amino acids as dictionary

    """
    properties_df = pd.read_csv(file_path, index_col='1-letter code')
    properties_series = properties_df['hydropathy index (Kyte-Doolittle method)']
    return dict(properties_series)


def calculate_hydropathy(sequence, amino_acid_properties):
    """Calculates the hydropathy of a amino acid sequence.

    Args:
        sequence: Amino acid sequence in 1-letter code as string.
        amino_acid_properties: Hydropathy of amino acids as dictionary.

    Returns: hydropathy as list

    """
    hydropathy = []
    for amino_acid in sequence:
        hydropathy.append(amino_acid_properties[amino_acid])
    return hydropathy


def rolling_mean(lst, window_size):
    """Calculates the rolling mean in a window.

    Args:
        lst: List containing the data for the rolling mean calculation.
        window_size: Size of the window as int.

    Returns: smoothed list

    """
    smoothed_lst = []
    window = deque(lst[:window_size - 1], maxlen=window_size)
    for value in lst[window_size - 1:]:
        window.append(value)
        window_sum = sum(window)
        smoothed_lst.append(window_sum / window_size)
    return smoothed_lst


def make_scatter_plot(hydropathy):
    """Generates a scatter plot with the hydropathy on y-axis and the position in sequence on x-axis.

    Args:
        hydropathy: List containing the hydropathy for each amino acid in the amino acid sequence.

    Returns: scatter plot

    """
    data = [
        go.Scatter(x=list(range(0, len(hydropathy))),
                   y=hydropathy
                   )
    ]
    layout = {
        'title': {
            'text': 'Human G-protein coupled receptor 183'
        },
        'yaxis': {
            'title': 'Hydropathy'
        },
        'xaxis': {
            'title': 'Postion'
        }
    }

    fig = go.Figure(data=data, layout=layout)
    fig.show()


if __name__ == '__main__':
    path_to_seq = '../data/G-protein_coupled_receptor.fasta'
    path_to_prop = '../data/amino_acid_properties.csv'
    aas = get_amino_acid_sequence(path_to_seq)
    aap = get_amino_acid_properties(path_to_prop)

    hydropathy_0 = calculate_hydropathy(aas, aap)
    make_scatter_plot(hydropathy_0)

    hydropathy_5 = rolling_mean(hydropathy_0, 5)
    make_scatter_plot(hydropathy_5)
