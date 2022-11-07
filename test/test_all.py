import scipy

from src.io.reader import read_csv
from omics.omics_integration import OmicsIntegration
from src.model.COBRAmodel import *
from src.experimental.ExpMatrix import *
from os.path import join
from graphics.plot import *
from stats.stats import *


def read_model(data_directory):
    model = MyModel(join(data_directory, "model.xml"), "e_Biomass__cytop")
    model.add_medium(join(data_directory, "media.xlsx"), "base_medium")
    model.reactions.e_Biomass_ht__cytop.bounds = (0, 0)
    model.reactions.R00019__chlo.bounds = (0, 0)
    model.reactions.R00019__mito.bounds = (0, 0)
    model.set_prism_reaction("PRISM_white_LED__extr")
    model.reactions.ATPm__cytop.bounds = (2.85, 2.85)
    return model


def experimental_data_processing(data_directory, filename, model):
    matrix = ExpMatrix(join(data_directory, filename), model)
    matrix.conditions = "Resume"
    matrix.remove_trials(["Resume", "area", "19", "21", "R21", "25", "PC1"])
    matrix.set_exponential_phases({"1": (2, 8), "2": (2, 8), "3": (2, 16), "4": (2, 10), "5": (2, 8), "6": (2, 8), "7": (2, 16), "8": (2, 16), "9": (2, 8), "10": (2, 8), "11": (2, 12),
                                   "12": (2, 10), "13": (2, 8) ,"14": (2, 8), "15": (2, 14), "16": (2, 14), "17": (2, 12), "18": (2, 12), "20": (2, 14),
                                   "22": (2, 14), "23": (2, 10), "24": (2, 10), "PC1": (2, 10), "PC2": (2, 14), "PC3": (2, 10), "PC4": (2, 10), "RPC1": (2, 10), "RPC2": (2, 10),
                                   "RPC3": (2, 10)})
    matrix.get_experimental_data(parameter='all')
    matrix.get_substrate_uptake_from_biomass("C", "CO2", header="C00011")
    matrix.get_substrate_uptake_from_biomass("P", "HPO4")
    matrix.get_substrate_uptake("[P] mmol", header='C00009')
    matrix.get_substrate_uptake("[N] mmol", header='C00244')
    matrix.save()
    return matrix


def simulations(matrix, model):
    complete_results, values_for_plot_carbon_limited, _ = simulation_for_conditions(model, matrix.conditions[["C00011"]], matrix.conditions[["growth_rate"]], save_in_file=True, filename="carbon_limited")
    _, values_for_plot_carbon_and_p_limited, _ = simulation_for_conditions(model, matrix.conditions[["C00011", "C00009"]], matrix.conditions[["growth_rate"]], save_in_file=True, filename="carbon_and_p_limited")
    _, values_for_plot_carbon_and_p_and_n_limited, _ = simulation_for_conditions(model, matrix.conditions[["C00011", "C00009", "C00244"]], matrix.conditions[["growth_rate"]], save_in_file=True, filename="carbon_and_p_and_n_limited")
    res = []
    for index, element in enumerate(values_for_plot_carbon_limited):
        new_value = values_for_plot_carbon_and_p_limited[index][1]
        new_value2 = values_for_plot_carbon_and_p_and_n_limited[index][1]
        res.append((element[0], element[1], new_value))
    df = pd.DataFrame(data=res, index=complete_results.keys(), columns=["Exp", "in silico (C limited)", "in silico (C & P limited)"])
    df.to_excel("simulation_summary.xlsx")
    barplot(df, to_show=False, path="simulation_summary.png")
    model.set_photoautotrophy()
    print(model.summary())
    model.set_heterotrophy()
    print(model.summary())
    model.set_mixotrophy()
    print(model.summary())


def stats(data_directory, filename):
    matrix = ExpMatrix(join(data_directory, filename))
    matrix.conditions = "Resume"
    matrix.conditions = matrix.conditions.rename({"[N] mmol": "N", "[P] mmol": "P", "Salinity g/L": "salinity", "Aeration rate": "aeration", 'growth_rate': 'umax', 'Productivity (g/L.h)': 'Pmax', 'Biomass (gDW/L)': 'biomass'}, axis=1)
    boxplot(matrix.conditions, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Pmax'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
                            , y_labels={'Pmax': 'Pmax $(g \cdot L^{-1} \cdot d^{-1})$'})
    stats = StatisticalAnalysis(matrix.conditions)
    anova_table, model = stats.manova('umax ~ P')
    anova_table, model = stats.anova('biomass ~ P')
    print(anova_table)
    hist(matrix.conditions, ['biomass'], title='Biomass', xlabel='Biomass $(g \cdot L^{-1})$', ylabel='Frequency')
    hist(matrix.conditions, ['Pmax'], title='Maximum Productivity', xlabel='Pmax $(g \cdot L^{-1} \cdot h^{-1})$', ylabel='Frequency')
    hist(matrix.conditions, ['umax'], title='Growth Rate', xlabel='Biomass $(h^{-1})$', ylabel='Frequency')
    qqplot(model, to_show=True)



def omics_integration(model):
    omics = OmicsIntegration('raw_counts.txt', samples_names={"PRJNA437866/SRR6825159/SRR6825159_Aligned.sortedByCoord.out.bam":"control_1",
                                                              "PRJNA437866/SRR6825160/SRR6825160_Aligned.sortedByCoord.out.bam": "control_2",
                                                               "PRJNA437866/SRR6825161/SRR6825161_Aligned.sortedByCoord.out.bam":"control_3",
                                                                "PRJNA437866/SRR6825162/SRR6825162_Aligned.sortedByCoord.out.bam":"nacl_1",
                                                                "PRJNA437866/SRR6825163/SRR6825163_Aligned.sortedByCoord.out.bam":"nacl_2",
                                                                "PRJNA437866/SRR6825164/SRR6825164_Aligned.sortedByCoord.out.bam":"nacl_3",
                                                                "PRJNA437866/SRR6825165/SRR6825165_Aligned.sortedByCoord.out.bam":"h2o2_1",
                                                                "PRJNA437866/SRR6825166/SRR6825166_Aligned.sortedByCoord.out.bam":"h2o2_2",
                                                                "PRJNA437866/SRR6825167/SRR6825167_Aligned.sortedByCoord.out.bam":"h2o2_3",
                                                                "PRJNA437866/SRR6825168/SRR6825168_Aligned.sortedByCoord.out.bam":"sorb_1",
                                                                "PRJNA437866/SRR6825169/SRR6825169_Aligned.sortedByCoord.out.bam":"sorb_2",
                                                                "PRJNA437866/SRR6825170/SRR6825170_Aligned.sortedByCoord.out.bam":"sorb_3",
                                                              }, model=model)

    # omics.get_getmm()
    # omics.get_degs()
    omics.getmm = read_csv("getmm.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics.degs = read_csv("degs.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    idx = omics.degs.index
    counts_degs = omics.getmm.loc[idx]
    # omics.counts.to_csv("counts.tsv", sep="\t")
    # omics.data.to_csv("data.tsv", sep="\t")
    # omics.sum_tech_reps()
    # omics.get_tpm()
    # omics.tpm.to_csv("tpm.tsv", sep="\t")
    g = clustermap(counts_degs)
    print()



    # omics.integrate(method='eFlux', parsimonious=True)
    # omics.integrate(method='GIMME', parsimonious=True, biomass="e_Biomass__cytop", growth_frac = 0.1)
    # omics.integrate(method='iMAT')
    # omics.save()


if __name__ == '__main__':
    data_directory = r"../data"
    model = read_model(data_directory)
    matrix = experimental_data_processing(data_directory, "Matriz- DCCR Dunaliella salina.xlsx", model)
    matrix = ExpMatrix(join(data_directory, "Matriz- DCCR Dunaliella salina_new.xlsx"))
    matrix.conditions = "Resume"
    simulations(matrix, model)
    biomass = Biomass("e_Biomass__cytop", "Biomass_exp_composition.xlsx")
    stats(data_directory, "Matriz- DCCR Dunaliella salina_new.xlsx")
    # omics_integration(None)