import os
import subprocess
from os.path import join
from pprint import pprint

import numpy as np
import pandas as pd
from bioinfokit.analys import norm
from mewpy.omics import ExpressionSet, eFlux, GIMME, iMAT

from src.bio.genes import Genes
from src.io.reader import read_csv
from src.io.writer import write_specific_models
from utils.utils import run


class OmicsIntegration:
    def __init__(self, filename: str, samples_names: dict = None, groups: list = None, model=None):
        self.model = model
        self.counts_file = join(os.getcwd(), filename)
        self.samples, self.samples_names = None, samples_names
        self._data, self._counts, self.tpm = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        self._genes = Genes()
        self.load()
        self.groups = groups or [e.split("_")[0] for e in self.samples]
        self.specific_models = {}

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, value):
        self._genes = Genes(batch=value)

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
        self.data = read_csv(self.counts_file, index_name='GeneID', index_col=0, comment='#', sep='\t')
        if self.samples_names:
            self.data = self.data.rename(columns=self.samples_names)
        if not self.samples:
            self.samples = sorted(list(set(self.data.columns) - {'Chr', 'Start', 'End', 'Strand', 'Length'}))
        self.counts = self.data[self.samples]
        self.genes = self.data['Length']

    def get_tpm(self):
        nm = norm()
        df_for_calculation = pd.concat([self.counts, self.genes.lengths], axis=1)
        nm.tpm(df=df_for_calculation, gl='Length')
        self.tpm = nm.tpm_norm


    def get_getmm(self, counts_file=None, data_path = None, output_file=None):
        if not output_file:
            output_file = join(os.getcwd(), "getmm.tsv")
        if not data_path:
            data_path = join(os.getcwd(), "data.tsv")
        if not os.path.exists(data_path):
            self.data.to_csv(data_path, sep='\t')
        if not counts_file:
            counts_file  = join(os.getcwd(), "counts.tsv")
        if not os.path.exists(output_file):
            self.counts.to_csv(counts_file, sep='\t')
        cmd = ["Rscript", join(os.getcwd(), "../src/omics/GeTMM.R"), counts_file , data_path, output_file, ','.join(self.groups)]
        run(cmd)
        self.getmm = read_csv(output_file, index_name='GeneID', index_col=0, comment='#', sep='\t')

    def get_degs(self, getmm_file=None, output_file=None, **kwargs):
        if not getmm_file:
            getmm_file = join(os.getcwd(), "getmm.tsv")
        if not output_file:
            output_file = join(os.getcwd(), "degs.tsv")
        cmd = ["Rscript", join(os.getcwd(), "../src/omics/DGE.R"), getmm_file, output_file, ','.join(self.groups)]
        run(cmd)
        self.degs = read_csv(output_file, index_name='GeneID', index_col=0, comment='#', sep='\t')




    def sum_tech_reps(self):
        for group in self.groups:
            idx = self.counts.columns.str.startswith(group)
            self.counts[group] = self.counts.iloc[:, idx].sum(axis=1)
        to_remove = [col for col in self.counts.columns if col not in self.groups]
        self.counts.drop(to_remove, axis=1, inplace=True)

    def integrate(self, method=None, samples=None, **kwargs):
        if not samples:
            samples = self.groups
        if method not in self.specific_models.keys():
            self.specific_models[method] = {}
        for sample in samples:
            expression = self.counts[sample].to_numpy()[:, np.newaxis]
            set_expression = ExpressionSet(self.genes.genes_ids, [sample], expression)
            method_callable = self.get_method(method)
            try:
                result = method_callable(self.model, expr=set_expression, condition=sample, **kwargs)
                self.specific_models[method][sample] = result
                pprint(result.dataframe.loc[result.dataframe.index=='EX_e_Biomass__dra'])
            except Exception as e:
                print(e)
                print(f"Error in {method} for {sample}")

    @staticmethod
    def get_method(method) -> callable:
        if method == 'eFlux':
            return eFlux
        elif method == 'GIMME':
            return GIMME
        elif method == 'iMAT':
            return iMAT

    def save(self):
        write_specific_models(self.specific_models, "context_specific_models.xlsx")
