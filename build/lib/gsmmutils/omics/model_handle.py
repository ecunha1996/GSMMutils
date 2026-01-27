import json
from os.path import abspath, join, dirname
import cobra
import os
import pandas as pd
from cobra.flux_analysis import find_blocked_reactions
from cobra.io import write_sbml_model
from cobra.manipulation.delete import prune_unused_metabolites


def print_model_details(cobra_model):
    """
    Function to print the details of the currently loaded COBRA model.

    Parameters
    ----------
    cobra_model : cobra.Model

    """
    transporters = []

    for reaction in cobra_model.reactions:
        if len(reaction.compartments) == 2:
            transporters.append(reaction.id)

    print('Total Reactions:', len(cobra_model.reactions))
    print('Reactions:', (len(cobra_model.reactions)) - len(transporters) - len(cobra_model.exchanges))
    print('Transporters:', len(transporters))
    print('Exchanges:', len(cobra_model.exchanges))


def load_model(model_path: str, consistent_model_path: str, medium_conditions) -> cobra.Model:
    """
    This function is used to load the model.

    Parameters
    ----------
    model_path : str
        The path to the model.
    consistent_model_path : str
        The path to the model without blocked reactions.

    Returns
    -------
    model : cobra.Model
        The loaded model.
    """
    if consistent_model_path and os.path.exists(consistent_model_path):
        model = cobra.io.read_sbml_model(consistent_model_path)
    else:
        model = cobra.io.read_sbml_model(model_path)
    for reaction_id, bound in medium_conditions.items():
        if reaction_id in model.reactions:
            model.reactions.get_by_id(reaction_id).bounds = bound
    # for demand in model.demands:
    #     if 'photon' not in demand.id:
    #         demand.bounds = (0, 0)
    if not consistent_model_path or not os.path.exists(consistent_model_path):
        model.remove_reactions(find_blocked_reactions(model, processes=30), remove_orphans=True)
        write_sbml_model(model, consistent_model_path)
    print_model_details(model)
    print(model.slim_optimize())
    return model


def sbml_model_reconstruction(original_model: cobra.Model, sample: str, integration_result_dict: dict, params, medium_conditions, model_results_path):
    """
    This function is used to reconstruct the model based on the integration results.

    Parameters
    ----------
    model_template: cobra.Model
        The COBRA model template.
    sample: str
        The sample name.
    integration_result_dict: dict
        The integration results.
    """
    #try:
    print(f"Reconstructing model for {sample}...")
    model_template = original_model.copy()
    model_template.objective = params['OBJECTIVE']
    reactions_to_deactivate = [reaction for reaction, value in
                               integration_result_dict[sample].items() if value is False and model_template.reactions.get_by_id(reaction).genes]
    print('Reactions to deactivate:', len(reactions_to_deactivate))
    for reaction in reactions_to_deactivate:
        model_template.remove_reactions([reaction], remove_orphans=True)

    for reaction_id, bound in medium_conditions[params['MEDIUM_NAME']].items():
        if reaction_id in model_template.reactions:
            model_template.reactions.get_by_id(reaction_id).bounds = bound
        else:
            print(reaction_id, 'exchange not found in the model.')
    model_template.remove_reactions(find_blocked_reactions(model_template, processes=min(os.cpu_count(), 30)), remove_orphans=True)
    model_template, _ = prune_unused_metabolites(model_template)
    model_name = sample + '.xml'
    print(model_template.slim_optimize())
    cobra.io.write_sbml_model(model_template, os.path.join(model_results_path, model_name))

    print(f'Model reconstruction for {sample} finished.')
    print_model_details(model_template)
