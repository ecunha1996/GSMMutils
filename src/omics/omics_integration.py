import warnings

import pandas as pd
from bioinfokit.analys import norm

from src.bio.genes import Genes
from src.io.reader import read_csv


class OmicsIntegration:
    def __init__(self, filename, samples=None, model=None):
        self.model = model
        self.counts_file = filename
        self.samples = samples
        self._data, self._counts = pd.DataFrame(), pd.DataFrame()
        self._genes = Genes()
        self.load()

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, value):
        to_load = []
        if self.model:
            genes_ids = [gene.id for gene in self.model.genes]
            for gene in value:
                if gene in genes_ids:
                    to_load.append(gene)
                else:
                    raise warnings.warn(f"Gene {gene} not found in model")
        else:
            to_load = value
        self._genes = Genes(batch=to_load)

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value

    @property
    def counts(self):
        return self._counts

    @counts.setter
    def counts(self, value):
        self._counts = value

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, value):
        self._samples = value

    def load(self):
        self.data = read_csv(self.counts_file, index_col=0, comment='#')
        if not self.samples:
            self.samples = list(set(self.data.columns) - {'Chr', 'Start', 'End', 'Strand', 'Length'})

        self.genes = self.data.index

    def tpm(self):
        nm = norm()
        df_for_calculation = pd.concat([self.counts, self.genes.lengths], axis=1)
        nm.tpm(df=df_for_calculation, gl='Length')
        tpm_df = nm.tpm_norm
        print(tpm_df.mean())
