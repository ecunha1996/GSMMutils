import itertools
import os
from os.path import join
from pprint import pprint

import numpy as np
import pandas as pd
# from bioinfokit.analys import norm
from cobra.flux_analysis import flux_variability_analysis as fva
from cobra.io import write_sbml_model
from mewpy.omics import ExpressionSet
from mewpy.omics import ExpressionSet, eFlux, GIMME, iMAT
from mewpy.simulation import get_simulator

from ..bio.genes import Genes
from ..graphics.plot import clustermap
from ..io import read_csv
from ..io import write_specific_models
from ..utils.utils import run, differential_reaction_capacity


class OmicsIntegration:
    def __init__(self, filename: str, samples_names: dict = None, groups: list = None, model=None):
        self.pathways_reaction_counts = None
        self.pathways_gene_counts = None
        self.degs = None
        self.getmm = None
        self.reaction_capacity = {}
        self.flux_change = {}
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
        self.data.index = self.data.index + "_1"
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

    def get_getmm(self, counts_file=None, data_path=None, output_file=None):
        if not output_file:
            output_file = join(os.getcwd(), "getmm.tsv")
        if not data_path:
            data_path = join(os.getcwd(), "data.tsv")
        if not os.path.exists(data_path):
            self.data.to_csv(data_path, sep='\t')
        if not counts_file:
            counts_file = join(os.getcwd(), "counts.tsv")
        if not os.path.exists(output_file):
            self.counts.to_csv(counts_file, sep='\t')
        cmd = ["Rscript", join(os.getcwd(), "../../src/gsmmutils/omics/GeTMM.R"), counts_file, data_path, output_file, ','.join(self.groups)]
        run(cmd)
        self.getmm = read_csv(output_file, index_name='GeneID', index_col=0, comment='#', sep='\t')

    def get_degs(self, getmm_file=None, output_file=None):
        if not getmm_file:
            getmm_file = join(os.getcwd(), "getmm.tsv")
        if not output_file:
            output_file = join(os.getcwd(), "degs.tsv")
        cmd = ["Rscript", join(os.getcwd(), "../../src/gsmmutils/omics/DGE.R"), getmm_file, output_file, ','.join(self.groups)]
        run(cmd)
        self.degs = read_csv(output_file, index_name='GeneID', index_col=0, comment='#', sep='\t')

    def set_pathways_counts_by_gene(self):
        pathways_map = self.model.get_genes_pathways_map()
        self.genes.set_pathways(pathways_map)
        counts_map = {}
        for gene in self.degs.index:
            if gene in pathways_map:
                for pathway in pathways_map[gene]:
                    if pathway not in counts_map:
                        counts_map[pathway] = 1
                    else:
                        counts_map[pathway] += 1
        self.pathways_gene_counts = pd.DataFrame.from_dict(counts_map, orient='index', columns=['Counts'])

    def set_pathways_counts_by_reaction(self, dataframe):
        counts_map = {}
        if not self.model.reactions_pathways_map: self.model.get_reactions_pathways_map()
        for reaction in dataframe.index:
            if reaction in self.model.reactions_pathways_map.keys():
                for pathway in self.model.reactions_pathways_map[reaction]:
                    if pathway not in counts_map:
                        counts_map[pathway] = 1
                    else:
                        counts_map[pathway] += 1
        self.pathways_reaction_counts = pd.DataFrame.from_dict(counts_map, orient='index', columns=['Counts'])

    def sum_tech_reps(self):
        for group in self.groups:
            idx = self.counts.columns.str.startswith(group)
            self.counts.loc[:, group] = self.counts.iloc[:, idx].sum(axis=1).values
        to_remove = [col for col in self.counts.columns if col not in self.groups]
        self.counts.drop(to_remove, axis=1, inplace=True)

    def integrate(self, method=None, samples=None, tool="mewpy", **kwargs):
        print(f"Integrating omics data with {method}")
        if not samples:
            samples = list(set(self.groups))
        if method not in self.specific_models.keys():
            self.specific_models[method] = {}
        if tool == "troppo":
            self.integrate_troppo(method, samples, **kwargs)
        else:
            for sample in samples:
                print("Integrating sample: {}".format(sample))
                expression = self.counts[sample].to_numpy()[:, np.newaxis]
                set_expression = ExpressionSet(self.genes.genes_ids, [sample], expression)
                method_callable = self.get_method(method)
                try:
                    if 'build_model' in kwargs and kwargs['build_model']:
                        res_sim, result = method_callable(self.model, expr=set_expression, condition=sample, **kwargs)
                        write_sbml_model(result.model, join(kwargs['data_path'], f"{sample}/Dsalina_{sample}_{method.lower()}.xml"))
                    else:
                        result = method_callable(self.model, expr=set_expression, condition=sample, **kwargs)
                    self.specific_models[method][sample] = result
                    if method == "eFlux":
                        constraints = {}
                        for reaction, bounds in result.simulation_constraints.items():
                            if bounds == (-1, 1) or bounds == (0, 1) or bounds == (-1, 0):
                                result.simulation_constraints[reaction] = self.model.reactions.get_by_id(reaction).bounds
                            constraints[reaction] = result.simulation_constraints[reaction]
                        as_df = pd.DataFrame.from_dict(constraints, orient='index', columns=['Lower', 'Upper'])
                        as_df.to_csv(join(os.getcwd(), "constraints_{}.tsv".format(sample)), sep='\t')
                        simulator = get_simulator(self.model, result.simulation_constraints)
                        sol = simulator.simulate()
                        sol.dataframe.to_csv(join(os.getcwd(), "results_{}_{}.tsv".format(method, sample)), sep='\t')
                        print(sol.status)
                        pprint(sol.dataframe.loc[result.dataframe.index == 'EX_e_Biomass__dra'])
                except Exception as e:
                    print(e)
                    print(f"Error in {method} for {sample}")
        print("#" * 100)

    @staticmethod
    def get_method(method) -> callable:
        if method == 'eFlux':
            return eFlux
        elif method == 'GIMME':
            return GIMME
        elif method == 'iMAT':
            return iMAT

    def save(self, filename="context_specific_models.xlsx"):
        write_specific_models(self.specific_models, filename)

    def get_clustermap_genes_by_pathway(self, omics_results=None):
        self.set_pathways_counts_by_gene()
        idx = self.degs.index
        count_degs = self.getmm.loc[idx]
        for pathway in self.pathways_gene_counts.index:
            if self.pathways_gene_counts.loc[pathway, 'Counts'] > 1:
                in_pathway = [gene_id for gene_id, value in self.model.genes_pathways_map.items() if pathway in value]
                clustermap(count_degs.loc[count_degs.index.isin(in_pathway)], title=pathway, to_show=False, path=f"./omics/clustermaps_degs/{pathway.replace(' ', '_').replace('/', '_')}.png")

    def get_clustermap_reactions_by_pathway(self, omics_results=None):
        for key, value in omics_results.items():
            value.columns = [key]
        new_df = pd.concat([value for key, value in omics_results.items()], axis=1)
        new_df = new_df.loc[(round(new_df.T, 7) != 0).any()]
        new_df = new_df.drop(["GIMME_h2o2", "GIMME_sorb"], axis=1)
        var = new_df.var(axis=1)
        variance = new_df.loc[var > 1e-5]
        self.set_pathways_counts_by_reaction(variance)
        for pathway in self.pathways_reaction_counts.index:
            if self.pathways_reaction_counts.loc[pathway, 'Counts'] > 1:
                in_pathway = [reaction_id for reaction_id, value in self.model.reactions_pathways_map.items() if pathway in value]
                clustermap(variance.loc[variance.index.isin(in_pathway)], title=pathway, to_show=False, path=f"./omics/\clustermaps_reactions/clustermaps_gimme/{pathway.replace(' ', '_').replace('/', '_')}.png")

    def integrate_troppo(self, method, samples, **kwargs):
        from ..omics.troppo import integration_pipeline
        for sample in samples:
            if "control" in sample:
                df = pd.DataFrame(self.counts[sample]).T
                median = 0
                integration_pipeline(dataset=df, dataset_name=sample, algorithm=method, threshold=median,
                                     thread_number=8, model=self.model.model.copy())

    def integrate_with_mewpy(self):
        pass

    def get_flux_change(self, combine_all=False, method_1=None, condition_1=None, method_2=None, condition_2=None, threshold: float = 0.1):
        from ..utils.utils import flux_change
        if combine_all:
            for method in self.specific_models.keys():
                for condition_1, condition_2 in itertools.product(self.specific_models[method].keys(), self.specific_models[method].keys()):
                    if condition_1 != condition_2:
                        constraints_1, constraints_2 = {}, {}
                        for reaction, flux in self.specific_models[method][condition_1].dataframe.to_dict(orient='index').items():
                            constraints_1[reaction] = flux['Flux rate']
                        for reaction, flux in self.specific_models[method][condition_2].dataframe.to_dict(orient='index').items():
                            constraints_2[reaction] = flux['Flux rate']
                        self.flux_change[f"{method}_{condition_1}_{condition_2}"] = flux_change(constraints_1, constraints_2, threshold=threshold)
        else:
            constraints_1, constraints_2 = {}, {}
            for reaction, flux in self.specific_models[method_1][condition_1].dataframe.to_dict(orient='index').items():
                constraints_1[reaction] = flux['Flux rate']
            for reaction, flux in self.specific_models[method_2][condition_2].dataframe.to_dict(orient='index').items():
                constraints_2[reaction] = flux['Flux rate']
            self.flux_change[f"{method_1}_{condition_1}_{method_2}_{condition_2}"] = flux_change(constraints_1, constraints_2, threshold)

    def get_reaction_capacity(self, condition, fva_solution):
        from ..utils.utils import reaction_capacity
        if not hasattr(self, "reaction_capacity"):
            self.reaction_capacity = {}
        self.reaction_capacity[condition] = reaction_capacity(fva_solution)

    def get_differential_reaction_capacity(self, method, condition_1, condition_2):
        fva_sol_1 = fva(self.specific_models[method][condition_1].model)
        self.get_reaction_capacity(condition_1, fva_sol_1)
        fva_sol_2 = fva(self.specific_models[method][condition_2].model)
        self.get_reaction_capacity(condition_2, fva_sol_2)
        return differential_reaction_capacity(self.reaction_capacity[condition_1], self.reaction_capacity[condition_2])
