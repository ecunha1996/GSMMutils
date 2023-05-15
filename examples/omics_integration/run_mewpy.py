import math
import os
from os.path import join

import numpy as np

from ExpGSMM.io import read_csv
from ExpGSMM.model.COBRAmodel import MyModel
from ExpGSMM.omics.omics_integration import OmicsIntegration
import pickle

def omics_integration(model):

    omics = OmicsIntegration('raw_counts.txt', samples_names={"SRR6825159_Aligned.sortedByCoord.out.bam":"control_1",
                                                              "SRR6825160_Aligned.sortedByCoord.out.bam": "control_2",
                                                               "SRR6825161_Aligned.sortedByCoord.out.bam":"control_3",
                                                                "SRR6825162_Aligned.sortedByCoord.out.bam":"nacl_1",
                                                                "SRR6825163_Aligned.sortedByCoord.out.bam":"nacl_2",
                                                                "SRR6825164_Aligned.sortedByCoord.out.bam":"nacl_3",
                                                                "SRR6825165_Aligned.sortedByCoord.out.bam":"h2o2_1",
                                                                "SRR6825166_Aligned.sortedByCoord.out.bam":"h2o2_2",
                                                                "SRR6825167_Aligned.sortedByCoord.out.bam":"h2o2_3",
                                                                "SRR6825168_Aligned.sortedByCoord.out.bam":"sorb_1",
                                                                "SRR6825169_Aligned.sortedByCoord.out.bam":"sorb_2",
                                                                "SRR6825170_Aligned.sortedByCoord.out.bam":"sorb_3",
                                                              }, model=model)

    # omics.get_getmm()
    # omics.get_degs()
    omics.getmm = read_csv("getmm.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    # omics.degs = read_csv("degs.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    # counts_degs = omics.counts.loc[omics.degs.index]
    # g = clustermap(counts_degs)
    omics.sum_tech_reps()
    omics.counts = omics.counts.applymap(lambda x: math.log2(x + 1))
    omics.counts.to_csv("counts_log2.tsv", sep="\t")
    import matplotlib.pyplot as plt
    omics.counts.plot.density()
    plt.xticks(np.arange(-10, 30, 1))
    plt.savefig(r"counts_density_nacl.png")
    tmp = omics.counts.loc[~(omics.counts == 0).all(axis=1)]
    print(tmp.describe())
    # omics.integrate(method='GIMME', biomass="e_Biomass__cytop", growth_frac=0.6, parsimonious=True)
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="nacl", method_2="GIMME")
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="h2o2", method_2="GIMME")
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="sorb", method_2="GIMME")
    # omics.integrate(method='eFlux', constraints= {"e_Biomass__cytop": (0.0015, 1000)})
    # omics.integrate(method='iMAT', constraints= {"e_Biomass__cytop": (0.15, 1000)}, cutoff=(15, 75))
    # omics.save("nacl_h2o2_sorb.xlsx")
    # pickle.dump(omics, open("nacl_h2o2_sorb/nacl_h2o2_sorb.pkl", "wb"))


def omics_integration_light_conditions(model):
    omics = OmicsIntegration('output.txt', samples_names={"SRR7984026Aligned.out.sam":"LL_1",
                                                                "SRR7984027Aligned.out.sam": "LL_2",
                                                                "SRR7984028Aligned.out.sam":"LL_3",
                                                                "SRR7984029Aligned.out.sam":"ML_1",
                                                                "SRR7984030Aligned.out.sam":"ML_2",
                                                                "SRR7984031Aligned.out.sam":"ML_3",
                                                                "SRR7984032Aligned.out.sam":"HL_1",
                                                                "SRR7984033Aligned.out.sam":"HL_2",
                                                                "SRR7984034Aligned.out.sam":"HL_3",
                                                              }, model=model)

    # omics.counts.to_csv(join(os.getcwd(),"counts_light.tsv"), sep="\t")
    # omics.data.to_csv(join(os.getcwd(),"data_light.tsv"), sep="\t")
    # omics.get_getmm(counts_file = join(os.getcwd(),"counts_light.tsv"), data_path = join(os.getcwd(),"data_light.tsv"), output_file=join(os.getcwd(),"getmm_light.tsv"))
    # omics.get_degs(getmm_file = join(os.getcwd(),"getmm_light.tsv"), output_file= join(os.getcwd(),"degs_light.tsv"))
    omics.getmm = read_csv("getmm_light.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics.degs = read_csv("degs_light.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics_degs = omics.counts.loc[omics.degs.index]
    # g = clustermap(omics_degs, to_show=False, path=join(os.getcwd(), "clustermap_light.png"))
    omics.sum_tech_reps()
    # omics.counts.index = omics.counts.index + "_1"
    # omics.counts.to_csv("counts_getmm_light.tsv", sep="\t")
    omics.counts = omics.counts.applymap(lambda x: math.log2(x + 1))
    # import matplotlib.pyplot as plt
    # omics.counts.plot.density()
    # plt.savefig(r"counts_density_light.png")
    omics.integrate(method='GIMME', biomass="e_Biomass__cytop", growth_frac=0.6, parsimonious=True)
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="nacl", method_2="GIMME")
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="h2o2", method_2="GIMME")
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="sorb", method_2="GIMME")
    # omics.model.reactions.e_Biomass__cytop.lower_bound = omics.model.optimize().objective_value * 0.1
    # omics.integrate(method='eFlux', max_exp=0.01)
    # omics.integrate(method='iMAT', constraints= {"e_Biomass__cytop": (0.15, 1000)}, cutoff=(15, 75))
    omics.save(filename="light/omics_light.xlsx")
    pickle.dump(omics, open("light/omics_light.pkl", "wb"))


if __name__ == '__main__':
    data_directory = r"../../data"
    model = MyModel(join(join(data_directory, 'models'), 'model_with_trials.xml'), "e_Biomass__cytop")
    model.add_medium(join(data_directory, "media.xlsx"), 'media_with_starch')
    os.chdir(rf"{data_directory}/omics")
    omics_integration(model)
    # omics_integration_light_conditions(model)