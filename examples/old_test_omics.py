# from troppo.methods_wrappers import GapfillWrapper
import math
import os

from GSMMutils.model.COBRAmodel import MyModel
from GSMMutils.io import read_csv
from GSMMutils.omics.omics_integration import OmicsIntegration
from GSMMutils.utils.utils import flux_change
from examples.test_all import read_model
import pandas as pd

def omics_integration(model):
    os.chdir(r"omics")
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

    omics.counts.to_csv("counts.tsv", sep="\t")
    omics.data.to_csv("data.tsv", sep="\t")
    omics.get_getmm()
    omics.get_degs()
    omics.getmm = read_csv("getmm.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    # omics.degs = read_csv("degs.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    # omics.counts.to_csv("counts.tsv", sep="\t")
    # g = clustermap(counts_degs)
    omics.sum_tech_reps()
    # omics.counts.index = omics.counts.index + "_1"
    omics.counts.to_csv("counts_getmm.tsv", sep="\t")
    omics.counts = omics.counts.applymap(lambda x: math.log2(x + 1))
    omics.integrate(method='GIMME', biomass="e_Biomass__cytop", growth_frac=0.7, parsimonious=True)
    omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="nacl", method_2="GIMME")
    omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="h2o2", method_2="GIMME")
    omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="sorb", method_2="GIMME")
    # omics.model.reactions.e_Biomass__cytop.lower_bound = omics.model.optimize().objective_value * 0.7
    # omics.integrate(method='eFlux')
    # omics.integrate(method='iMAT')
    # omics.save()

    model.get_reactions_pathways_map()

    pathway_counts_nacl_gimme = count_reaction_by_pathway(omics.flux_change['GIMME_control_GIMME_nacl'], model)
    pathway_counts_sorb_gimme = count_reaction_by_pathway(omics.flux_change['GIMME_control_GIMME_h2o2'], model)
    pathway_counts_h2o2_gimme = count_reaction_by_pathway(omics.flux_change['GIMME_control_GIMME_sorb'], model)

    print(pathway_counts_nacl_gimme)
    print("#" * 300)
    print(pathway_counts_h2o2_gimme)
    print("#" * 300)
    print(pathway_counts_sorb_gimme)

    # omics_results = read_matrix("omics/context_specific_models.xlsx", index_col=0, sheet_name=None)
    # omics.get_clustermap_genes_by_pathway()
    # omics.get_clustermap_reactions_by_pathway(omics_results)


def omics_integration_light_conditions(model):
    os.chdir(r"omics")
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
    import matplotlib.pyplot as plt
    omics.counts.plot.density()
    plt.savefig(r"counts_density_light.png")
    omics.model.reactions.DM_C00369__chlo.bounds = (-0.001, 1000)
    omics.integrate(method='GIMME', biomass="e_Biomass__cytop", growth_frac=0.6, parsimonious=True)
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="nacl", method_2="GIMME")
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="h2o2", method_2="GIMME")
    # omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="sorb", method_2="GIMME")
    # omics.model.reactions.e_Biomass__cytop.lower_bound = omics.model.optimize().objective_value * 0.1
    # omics.integrate(method='eFlux', max_exp=0.01)
    # omics.integrate(method='iMAT')
    omics.save(filename="light/omics_light.xlsx")

    # model.get_reactions_pathways_map()
    #
    # pathway_counts_nacl_gimme = count_reaction_by_pathway(omics.flux_change['GIMME_control_GIMME_nacl'], model)
    # pathway_counts_sorb_gimme = count_reaction_by_pathway(omics.flux_change['GIMME_control_GIMME_h2o2'], model)
    # pathway_counts_h2o2_gimme = count_reaction_by_pathway(omics.flux_change['GIMME_control_GIMME_sorb'], model)
    #
    # print(pathway_counts_nacl_gimme)
    # print("#" * 300)
    # print(pathway_counts_h2o2_gimme)
    # print("#" * 300)
    # print(pathway_counts_sorb_gimme)

    # omics_results = read_matrix("omics/context_specific_models.xlsx", index_col=0, sheet_name=None)
    # omics.get_clustermap_genes_by_pathway()
    # omics.get_clustermap_reactions_by_pathway(omics_results)


def get_flux_changes(model: MyModel):
    simulation_data = pd.read_excel(r"omics/context_specific_models.xlsx", sheet_name=None)
    fluxes_as_dict = {}
    for key, value in simulation_data.items():
        value.index = value["Reaction ID"]
        value = value.drop(["Reaction ID"], axis=1)
        temp_dict = value.to_dict(orient='index')
        for temp_key, temp_value in temp_dict.items():
            temp_dict[temp_key] = temp_value["Flux rate"]
        fluxes_as_dict[key] = temp_dict

    gimme_control_nacl = flux_change(fluxes_as_dict["GIMME_control"], fluxes_as_dict["GIMME_nacl"])
    gimme_control_h2o2 = flux_change(fluxes_as_dict["GIMME_control"], fluxes_as_dict["GIMME_h2o2"])
    gimme_control_sorb = flux_change(fluxes_as_dict["GIMME_control"], fluxes_as_dict["GIMME_sorb"])

    # enzymes_count = {}
    # for key, value in gimme_control_nacl.items():
    #     reaction = model.reactions.get_by_id(key)
    #     if 'ec-code' in reaction.annotation.keys():
    #         if type(reaction.annotation['ec-code']) == str: reaction.annotation['ec-code'] = [reaction.annotation['ec-code']]
    #         for ec in reaction.annotation['ec-code']:
    #             if ec not in enzymes_count.keys():
    #                 enzymes_count[ec] = 1
    #             else:
    #                 enzymes_count[ec] += 1
    # print(enzymes_count)
    model.get_reactions_pathways_map()

    pathway_counts_nacl_gimme = count_reaction_by_pathway(gimme_control_nacl, model)
    pathway_counts_h2o2_gimme = count_reaction_by_pathway(gimme_control_h2o2, model)
    pathway_counts_sorb_gimme = count_reaction_by_pathway(gimme_control_sorb, model)

    print(pathway_counts_nacl_gimme)
    print("#" * 300)
    print(pathway_counts_h2o2_gimme)
    print("#" * 300)
    print(pathway_counts_sorb_gimme)


def count_reaction_by_pathway(flux_changes: dict, model: MyModel):
    pathway_counts = {}
    for reaction in flux_changes.keys():
        for pathway in model.reactions_pathways_map[reaction]:
            if pathway in pathway_counts.keys():
                pathway_counts[pathway] += 1
            else:
                pathway_counts[pathway] = 1
    return pathway_counts

def test_csm():
    model = MyModel(r"C:\Users\Bisbii\PythonProjects\ExpAlgae\data\omics\control_tinit_t1.xml", 'e_Biomass__cytop')
    print(model.objective.expression)
    print(model.optimize())
    model.biomass_reaction = "e_Biomass__cytop"
    print(model.exchanges.EX_C00001__dra.bounds)
    print(model.reactions.T_H2O__plas.bounds)
    print(model.test_reaction("EX_C00001__dra"))

def sample():
    model = MyModel(r"C:\Users\Bisbii\PythonProjects\ExpAlgae\data\model_with_trials.xml", "e_Biomass__cytop")
    for reaction in model.reactions:
        if reaction.lower_bound < -1000:
            reaction.lower_bound = -1000
        if reaction.upper_bound > 1000:
            reaction.upper_bound = 1000
    #
    # print(model.objective.expression)
    reactions_to_remove = pd.read_csv(r"C:\Users\Bisbii\PythonProjects\ExpAlgae\data\omics\tinit_t18.0.csv", index_col=0)
    reaction_ids = []
    for reaction in reactions_to_remove.columns:
        res = reactions_to_remove[reaction].values[0]
        if res and model.reactions.get_by_id(reaction).genes and "EX" not in reaction:
            reaction_ids.append(reaction)
    print(len(reaction_ids))
    for reaction in reaction_ids:
        copy = model.copy()
        # model.reactions.get_by_id(reaction).knock_out()
        model.remove_reactions([reaction])
        sol = model.optimize()
        if sol.objective_value < 0.1 or sol.status != "optimal":
            model = copy
    print(model.optimize().objective_value)
    print(len(model.reactions))
    # model.write(r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\control_tinit_t18.1.xml")
    # cobra_model = cobra.io.read_sbml_model(r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\control_tinit_t18.1.xml")
    # print(len(cobra_model.reactions))
    # print(cobra_model.optimize().objective_value)
    # res = model.sample(constraints={"e_Biomass__cytop": (0.1, 1000)})
    # print(res)

if __name__ == '__main__':
    data_directory = r"../data"
    model = read_model(data_directory, 'model_with_trials.xml')
    omics_integration(model)
    # omics_integration_light_conditions(model)
    # get_flux_changes(model)
    # sample()
    # test_csm()
    # data = pd.read_csv("omics/tinit_t1.csv", index_col=0)
    # as_dict = data.to_dict(orient = 'index')
    # new_dict = as_dict['control']
    # with model as temp_model:
    #     temp_model.objective = 'e_Biomass__cytop'
    #
    #     reactions_to_deactivate = [reaction for reaction, value in
    #                                new_dict.items() if not value and temp_model.reactions.get_by_id(reaction).genes and "EX" not in reaction]
    #
    #     for reaction in reactions_to_deactivate:
    #         model.reactions.get_by_id(reaction).knock_out()
    #     non_consumed = {b.id for b in temp_model.boundary if b.bounds[0] >= 0}
    #     ls_override = {'produced': {"e_Biomass__cytop"},
    #                    'non_consumed': [list(temp_model.reactions.get_by_id(r).metabolites)[0].id for r in non_consumed]}
    #     gpfl_wrapper = GapfillWrapper(model=temp_model.model)
    #     sol = set(gpfl_wrapper.run(avbl_fluxes=reactions_to_deactivate, algorithm='efm')[0])
    #     print(sol)
    #     print(len(reactions_to_deactivate))
    #     reactions_to_deactivate = set(reactions_to_deactivate) - sol
    #     print(len(reactions_to_deactivate))
    #     for reaction in reactions_to_deactivate:
    #         copy = model.copy()
    #         copy.remove_reactions([reaction], remove_orphans=True)
    #         if copy.optimize().objective_value < 0.1:
    #             print(temp_model.reactions.get_by_id(reaction).id)
    #     print(temp_model.summary())
