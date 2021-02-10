import requests
import pandas as pd
import plotly.graph_objects as go
from collections import deque


def make_lookup():
    """Generates a dictionary with property-values of each amino acid for different physical properties.

    Returns: a dictionary with property-values of each amino acid for different physical properties

    """
    names_of_properties = [
        'Molecular Weight', 'Residue Weight', 'pka1', 'pka2', 'pkaX', 'pI',
        'hydropathy index (Kyte-Doolittle method)', 'Accessible surface'
    ]
    df_properties = pd.read_csv('../data/amino_acid_properties.csv', index_col='1-letter code')
    lookup = {}
    for prop_name in names_of_properties:
        property_series = df_properties[prop_name]
        lookup[prop_name] = dict(property_series)
    return lookup


class Protein:
    def __init__(self, lookup):
        self.lookup = lookup
        self.Id = self.url = self.fasta = self.sequence = self.prop_values = self.prop = self.window_size = None

    def get_data(self, Id):
        """Extracts the amino acid sequence for a protein Id.

        Args:
            Id: Protein Id as string.

        Returns: amino acid sequence as list

        """
        self.Id = Id
        self.url = 'https://www.uniprot.org/uniprot/' + self.Id + '.fasta'
        res = requests.get(self.url)
        self.fasta = res.text
        self.sequence = ''
        for line in self.fasta.split('\n'):
            if not line.startswith('>'):
                self.sequence += line.strip()
        return self.sequence

    def map_property(self, selected_property='hydropathy index (Kyte-Doolittle method)'):
        """Gives the property values for the input amino acid sequence.

        Args:
            selected_property: Physical property for which the property values are calculated.

        Returns: list of property values

        """
        self.prop = selected_property
        self.window_size = None
        prop_dic = self.lookup[self.prop]
        self.prop_values = [prop_dic[amino_acid] for amino_acid in self.sequence]
        return self.prop_values

    def rolling_mean(self, selected_property='hydropathy index (Kyte-Doolittle method)', window_size=5):
        """Gives the property values calculated by rolling mean for the input amino acid sequence.

        Args:
            selected_property: Physical property for which the property values are calculated.
            window_size: Size of the window as int.

        Returns: List of the property values calculated by rolling mean for the input amino acid sequence

        """
        prop_values_unsmoothed = self.map_property(selected_property)
        self.prop_values = []
        self.window_size = window_size

        window = deque(prop_values_unsmoothed[:window_size - 1], maxlen=window_size)
        for value in prop_values_unsmoothed[window_size - 1:]:
            window.append(value)
            window_sum = sum(window)
            self.prop_values.append(window_sum / window_size)
        return self.prop_values

    def make_plot(self):
        """Generates a scatter plot for the property values.

        Returns: scatter plot for the property values

        """
        data = [
            go.Scatter(x=list(range(0, len(self.prop_values))),
                       y=self.prop_values
                       )
        ]
        plot_title = self.Id
        if self.window_size is not None:
            plot_title += f' with window size of {self.window_size}'
        layout = {
            'title': {
                'text': plot_title
            },
            'yaxis': {
                'title': self.prop
            },
            'xaxis': {
                'title': 'Postion'
            }
        }

        fig = go.Figure(data=data, layout=layout)
        fig.show()


if __name__ == '__main__':
    ID = 'P32249'
    look_up = make_lookup()
    # prop = 'pI'
    prop = 'Accessible surface'
    # prop = 'hydropathy index (Kyte-Doolittle method)'

    prot = Protein(look_up)
    seq = prot.get_data(ID)
    prop_values = prot.map_property(prop)
    # prop_values_smooted = prot.rolling_mean(prop, 15)
    prot.make_plot()
