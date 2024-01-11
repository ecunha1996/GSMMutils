import traceback
import warnings
from os.path import join

from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.readers.generic import TabularReader

from gsmmutils.model.COBRAmodel import MyModel
from gsmmutils.utils.utils import enable_print

warnings.filterwarnings("ignore")

import cobra
import pandas as pd
import re
#
from cobamp.utilities.parallel import batch_run


#
#
# protected_reactions = ['Biomass_vvinif2021', 'e-Protein_vvinif2021', 'e-RNA_vvinif2021', 'e-DNA_vvinif2021',
#                        'e-Carbohydrates_vvinif2021', 'e-Lipids_vvinif2021', 'e-Cofactor_vvinif2021',
#                        'e-CellWall_vvinif2021', 'Maintenance_vvinif2021']
#
# properties_init = {'essential_reactions': protected_reactions}


def print_model_details(cobra_model):
    """
    Function to print the details of the currently loaded COBRA model.

    Parameters
    ----------
    cobra_model : cobra.Model

    """
    enable_print()

    transporters = []

    for reac in cobra_model.reactions:
        if len(reac.compartments) == 2:
            transporters.append(reac.id)

    print('Total Reactions:', len(cobra_model.reactions))
    print('Reactions:', (len(cobra_model.reactions)) - len(transporters) - len(cobra_model.exchanges))
    print('Transporters:', len(transporters))
    print('Exchanges:', len(cobra_model.exchanges))


#
def reconstruction_function(omics_container, parameters: dict):
    """
    This function is used to run the reconstruction algorithm.

    Parameters
    ----------
    omics_container : pandas.DataFrame
        The omics data set.
    parameters : dict
        The parameters to be used for the reconstruction algorithm.

    Returns
    ----------
    rec_wrapper : Reconstruction Wrapper object with model and omics data.
    """

    def integration_fx(data_map):
        """
        Custom integration function for the reconstruction algorithm.
        This selects reactions with scores higher than the threshold, while saving the biomass reaction.

        """

        return [[reaction for reaction, score in data_map.get_scores().items()
                 if (score is not None and score > threshold) or reaction in parameters['protected']]]

    def score_apply(data_map):
        dm = {k: 0 if v is None else (min(v, 10) - threshold) if k not in parameters['protected'] else 20
              for k, v in data_map.items()}
        return dm

    threshold, rec_wrapper, method = [parameters[parameter] for parameter in
                                      ['threshold', 'reconstruction_wrapper', 'algorithm']]

    # noinspection PyBroadException
    try:
        if method == 'fastcore':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method,
                                              integration_strategy=('custom', [integration_fx]), solver='CPLEX')

        elif method == 'tinit':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method,
                                              integration_strategy=('continuous', score_apply), solver='CPLEX')

    except:
        traceback.print_exc()

        return {r: False for r in rec_wrapper.model_reader.r_ids}


#
#
#
def troppo_integration(template_model, omics_dataset, algorithm: str, threshold: float,
                       thread_number: int, details: list, troppo_res_path: str):
    """
    This function is used to run the Troppo's integration algorithms.

    Parameters
    ----------
    template_model: cobra.Model
        The template model.
    omics_dataset: pandas.DataFrame
        The omics data set.
    algorithm: str
        The algorithm to be used.
    threshold: float
        The threshold to be used.
    thread_number: int
        The number of threads to be used for the omics integration.
    details: list
        Details of the integration performed.
        details[0] = model_name
        details[1] = algorithm
        details[2] = threshold
    troppo_res_path : str
        The path to save the results of the Troppo integration.

    Returns
    -------
    integration_results: dict
        Dataframe containing the results of the omics integration.
        Each sample as a dictionary containing a boolean value for each reaction.

    """

    patt = re.compile('__COBAMPGPRDOT__[0-9]')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    omics_data = TabularReader(path_or_df=omics_dataset,
                               omics_type='transcriptomics',
                               nomenclature='custom',
                               ).to_containers()
    reconstruction_wrapper = ReconstructionWrapper(template_model,
                                                   gpr_gene_parse_function=replace_alt_transcripts)

    enable_print()
    print('Reconstruction Wrapper Finished.')

    parameters = {'threshold': threshold, 'reconstruction_wrapper': reconstruction_wrapper, 'algorithm': algorithm, "flux_threshold": 1e-6,  #
                  "protected": ["e_Biomass__cytop"]}  # reaction.id for reaction in template_model.reactions
    # if not reaction.genes

    batch_fastcore_res = batch_run(reconstruction_function, omics_data, parameters, threads=thread_number)

    result_dict = dict(zip([sample.condition for sample in omics_data], batch_fastcore_res))

    integration_result = pd.DataFrame.from_dict(result_dict, orient='index')
    integration_result.to_csv(troppo_res_path)

    enable_print()
    print('Omics Integration with the %s method Finished.' % details[1])

    return result_dict


#
#
# def define_medium_conditions(model_template: cobra.Model):
#     """
#     Function to define the medium conditions for the context-specific models.
#
#     Parameters
#     ----------
#     model_template: cobra.Model
#         The template model.
#
#     Returns
#     ----------
#     model_template: cobra.Model
#         The template model with the defined medium conditions.
#
#     """
#
#     enable_print()
#
#     print('Defining medium conditions.')
#
#     medium_conditions = {}
#
#     for reaction_id, bound in medium_conditions.items():
#         if reaction_id in model_template.reactions:
#             model_template.reactions.get_by_id(reaction_id).bounds = bound
#         else:
#             print(reaction_id, 'exchange not found in the model.')
#
#     block_print()
#
#     return model_template
#
#
def reconstruct_context_specific_models(model_template: cobra.Model, integration_result_dict: dict, details: list, output_dir: str):
    """
    Function to obtain the context-specific models from the troppo integration results.

    Parameters
    ----------
    model_template: cobra.Model
        The template model.
    integration_result_dict: dict
        The dictionary containing the results of the omics integration.
        Each sample as a dictionary containing a boolean value for each reaction.
    details: list
        Details of the integration performed.
        details[0] = model_name
        details[1] = algorithm
        details[2] = threshold

    """
    as_df = pd.DataFrame.from_dict(integration_result_dict, orient='index')
    as_df.to_csv(join(output_dir, 'results.csv'))
    for sample in integration_result_dict.keys():
        enable_print()
        print('-----------------------------------------')
        print(sample)
        print('-----------------------------------------')
        with model_template as temp_model:
            temp_model.objective = 'e_Biomass__cytop'

            reactions_to_deactivate = [reaction for reaction, value in
                                       integration_result_dict[sample].items() if not value and temp_model.reactions.get_by_id(reaction).genes and "EX" not in reaction]

            to_remove = []
            gapfilling_reactions = []
            enable_print()
            print('Deactivating reactions.:' + str(len(reactions_to_deactivate)))
            print(temp_model.optimize())
            temp_model.remove_reactions(reactions_to_deactivate)
            # for r in reactions_to_deactivate:
            #     reaction = temp_model.reactions.get_by_id(r)
            #     copy_model = temp_model.copy()
            #     temp_model.remove_reactions([reaction])
            #     sol = temp_model.optimize()
            #     if sol.objective_value > 0.1:
            #         print(f'Reaction {r} is ok.')
            #     else:
            #         # gapfilling_reactions.append(reaction.id)
            #         temp_model = copy_model
            #     print(temp_model.optimize())
            # with open(f"omics/gapfilling_reactions_{sample}.txt", "w") as f:
            #     f.write("\n".join(gapfilling_reactions) + "\n")
            # with open(f"omics/reactions_to_remove_{sample}.txt", "w") as f:
            #     f.write("\n".join([r.id for r in to_remove]) + "\n")

            print(temp_model.optimize())
            print('Model for the %s sample has been generated.' % sample)

            print_model_details(temp_model)

            cobra.io.write_sbml_model(temp_model, join(output_dir, f'{sample}_{details[1]}_{details[2]}.xml'))


#
def integration_pipeline(dataset: pd.DataFrame, dataset_name: str, algorithm: str, threshold: float, thread_number: int, model: MyModel, output_dir="./"):
    """
    This function is used to run an integration pipeline with troppo.
    This function exports a context-specific SBML model for each sample of the dataset.

    Parameters
    ----------
    dataset: pd.DataFrame
        The name of the dataset.
    algorithm: str
        The algorithm to be used.
    threshold: float
        The threshold to be used.
    thread_number: int
        The number of threads to be used for the omics integration.

    """
    integration_details = [dataset_name, algorithm, threshold]

    directories = dict(Troppo_results=join(output_dir, f'{integration_details[1]}_t{integration_details[2]}.csv'))

    # if algorithm == 'fastcore':
    #     directories = {'Model': 'Marta/vvinif2021_v500_comparts_transports_noconstraints_noblocked.xml',
    #                    'Omics': 'Marta/GSE36128_final_transp.csv',
    #                    'Troppo_results': 'Marta/fastcore_results.csv'}
    #
    # elif algorithm == 'tinit':
    #     directories = {'Model': 'Marta/vvinif2021_v500_comparts_transports_noconstraints_noblocked.xml',
    #                    'Omics': 'Marta/GSE36128_final_transp.csv',
    #                    'Troppo_results': 'Marta/tinit_results.csv'}

    # define_medium_conditions(model)

    print_model_details(model)
    result = troppo_integration(model, dataset, algorithm, threshold, thread_number,
                                integration_details, directories['Troppo_results'])
    reconstruct_context_specific_models(model, result, integration_details, output_dir=output_dir)
