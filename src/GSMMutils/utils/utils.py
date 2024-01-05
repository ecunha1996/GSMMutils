import json
import math
import os
import subprocess
import sys
from os.path import join, abspath, dirname

import numpy as np
import pandas as pd
from cobra import Metabolite
from joblib import Parallel, delayed
from scipy.stats import linregress

from GSMMutils.experimental.BiomassComponent import BiomassComponent

CONFIG_PATH = abspath(join(dirname(__file__), '../../../config'))


def get_login_info(server: str):
    try:
        if os.path.exists(join(CONFIG_PATH, 'server_connection.json')):
            with open(join(CONFIG_PATH, 'server_connection.json')) as f:
                data = json.load(f)
                host = data[server]['host']
                username = data[server]['username']
                password = data[server]['password']
        else:
            host = input("Enter host: ")
            username = input("Enter username: ")
            password = input("Enter password: ")
        return host, username, password
    except Exception as e:
        print(e)
        sys.exit(1)


def get_productivity(data, time=48):
    previous_dw = 0
    productivity_g_l = [0]
    for i, row in data.iterrows():
        if int(i) == 0:
            previous_dw = row["DW"]
        if int(i) > 0:
            productivity_g_l.append((row["DW"] - previous_dw) / time)
            previous_dw = row["DW"]
    return productivity_g_l


def get_fixation(data, m_substrate, m_element, carbon_in_biomass):
    fixation = []
    for i, row in data.iterrows():
        fixation.append(row["Productivity (g/L.h)"] * carbon_in_biomass * (m_substrate / m_element))
    return fixation


def get_uptake(data, m_substrate, substrate):
    previous_dw = 0
    uptake = [0]
    for i, row in data.iterrows():
        if int(i) == 0:
            previous_dw = row["DW"]
        if int(i) > 0:
            avv = (row["DW"] + previous_dw) / 2
            uptake.append((row[f"{substrate} Fixation (g{substrate}/L.h)"] / avv) * 1000 / m_substrate)
            previous_dw = row["DW"]
    return uptake


def get_data_at_phase(data, exponential_phase):
    time_points = np.arange(exponential_phase[0] + 1, exponential_phase[1] + 1)
    sub_data = data.loc[data.index.astype(int).isin(time_points)]
    return sub_data


def get_average_uptake(data, exponential_phase, substrate) -> float:
    sub_data = get_data_at_phase(data, exponential_phase)
    return np.mean(sub_data[f"{substrate} uptake (mmol{substrate}/gDW.h)"]).round(4)


def get_element_in_biomass(model, element, biomass_reaction):
    res = get_biomass_mass(model, biomass_reaction)
    print(res)
    percentage_map = {}
    for key in res[1]:
        percentage_map[key] = round(res[1][key] * Metabolite(formula=key).formula_weight / 1000, 3)
    if round(sum(percentage_map.values()), 3) != 1: print(
        f"Error! Sum of Elemental mass percentage is different from 1! ({round(sum(percentage_map.values()), 3)})")
    return percentage_map[element]


def get_molecular_weight(formula):
    return Metabolite(formula=formula).formula_weight


def get_growth_rate_from_slope(data, exponential_phase):
    sub_data = get_data_at_phase(data, exponential_phase)
    a = sub_data.index.astype(int).to_list()
    b = np.log(sub_data["DW"].astype(float).to_list())
    return round(linregress(a, b).slope, 3)


def get_growth_rate(data, exponential_phase):
    a = np.log(data.loc[str(exponential_phase[0])]["DW"])
    b = np.log(data.loc[str(exponential_phase[1])]["DW"])
    return round((b - a) / (exponential_phase[1] - exponential_phase[0]), 3)


def get_maximum_productivity(data, exponential_phase):
    a = data.loc[str(exponential_phase[0])]["DW"]
    b = data.loc[str(exponential_phase[1])]["DW"]
    return round((b - a) / (exponential_phase[1] - exponential_phase[0]), 3)


def run(cmd):
    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.stdout.decode("utf-8"), process.stderr.decode("utf-8")
    if stdout:
        print(stdout)
    if stderr:
        print("\n\nError: ", stderr)


def get_precursors(macromolecule, temp_precursor, model):
    to_ignore = ["C00001__cytop"]
    precedent_reaction = [reaction for reaction in temp_precursor.reactions if
                          "Biomass" not in reaction.id and reaction.id.startswith("e_") and reaction.metabolites[
                              temp_precursor] > 0]
    if len(precedent_reaction) > 1:
        raise Exception(f"More than one reaction found for {temp_precursor.id}")
    if len(precedent_reaction) == 0:
        return temp_precursor.id
    else:
        if macromolecule.id not in to_ignore:
            for precursor in precedent_reaction[0].reactants:
                precursor = get_precursors(macromolecule, precursor, model)
                if type(precursor) == str:
                    if temp_precursor.id in model.biomass_components:
                        parent_component = model.biomass_components[temp_precursor.id]
                    else:
                        parent_component = BiomassComponent(temp_precursor, precedent_reaction[0].metabolites[
                            model.metabolites.get_by_id(precursor)], model.biomass_components[macromolecule.id])
                        model.biomass_components[temp_precursor.id] = parent_component
                    BiomassComponent(model.metabolites.get_by_id(precursor),
                                     precedent_reaction[0].metabolites[model.metabolites.get_by_id(precursor)],
                                     parent_component)


def update_st(stoichiometries, new_value):
    for key, value in new_value.items():
        stoichiometries[key] = value
    mysum = sum(stoichiometries.values())
    for key, value in stoichiometries.items():
        stoichiometries[key] = value / abs(mysum)
    return stoichiometries


def get_biomass_mass(model, biomass_reaction=None, lipid_subreactions = None):
    if lipid_subreactions is None:
        lipid_subreactions = [model.reactions.e_TAG__lip,
         model.reactions.e_DAG__er, model.reactions.e_DGTS__er, model.reactions.e_PE__cytop, model.reactions.e_PC__cytop,
         model.reactions.e_PI__er, model.reactions.e_PG__chlo, model.reactions.e_DGDG__chlo,
         model.reactions.e_SQDG__chlo, model.reactions.e_MGDG__chlo, model.reactions.e_CL__mito,
         model.reactions.e_FFA__cytop]
    def get_sum_of_reaction(reaction, stoichiometry, ignore_water, elementar_counter):
        counter = 0
        for reactant in reaction.reactants:
            copy = reactant.copy()
            copy.elements = {key: value for key, value in copy.elements.items() if key != "T"}
            counter += abs(reaction.metabolites[reactant]) * copy.formula_weight * stoichiometry
            for key in elementar_counter.keys():
                if key in reactant.elements:
                    elementar_counter[key] += abs(reaction.metabolites[reactant]) * reactant.elements[
                        key] * stoichiometry
        for product in reaction.products:
            if product.id != "C00001__cytop" or not ignore_water:
                copy = product.copy()
                copy.elements = {key: value for key, value in copy.elements.items() if key != "T"}
                counter -= reaction.metabolites[product] * copy.formula_weight * stoichiometry
                for key in elementar_counter.keys():
                    if key in product.elements:
                        elementar_counter[key] -= abs(reaction.metabolites[product]) * product.elements[
                            key] * stoichiometry
        return round(counter / 1000, 5), elementar_counter

    def parse_lipids(current_biomass_reaction, current_element_counter):
        current_counter = 0
        lipid_stoichiometry = abs(current_biomass_reaction.metabolites[model.metabolites.e_Lipid__cytop])
        lipids_reaction = model.reactions.e_Lipid__cytop

        for reaction in lipid_subreactions:
            for reactant in reaction.reactants:
                current_counter += abs(reaction.metabolites[reactant] * reactant.formula_weight * lipids_reaction.metabolites[
                    reaction.products[0]] * lipid_stoichiometry)
                for key in current_element_counter.keys():
                    if key in reactant.elements:
                        current_element_counter[key] += abs(
                            reaction.metabolites[reactant] * reactant.elements[key] * lipids_reaction.metabolites[
                                reaction.products[0]] * lipid_stoichiometry)
        return current_counter, current_element_counter

    if not biomass_reaction:
        biomass_reaction = model.bio_reaction
    elif type(biomass_reaction) == str:
        biomass_reaction = model.reactions.get_by_id(biomass_reaction)
    element_counter = {"C": 0, "N": 0, "O": 0, "P": 0, "S": 0, "H": 0}
    counter, element_counter = parse_lipids(biomass_reaction, element_counter)
    c = 0
    c += counter / 1000
    to_ignore = ["C00002__cytop", "C00001__cytop", "e_Lipid__cytop"]
    for reactant in biomass_reaction.reactants:
        if reactant.id not in to_ignore:
            for reaction in reactant.reactions:
                if "biomass" not in reaction.id.lower() and reaction.bounds != (0, 0):
                    go = True
                    if "trial" in reaction.id:
                        trial_number = reaction.id.split("trial")[1].split("_")[0]
                        if "trial" in biomass_reaction.id:
                            biomass_trial = biomass_reaction.id.split("trial")[1].split("_")[0]
                            if trial_number != biomass_trial:
                                go = False
                        else:
                            go = False
                    if go:
                        print(reaction.id)
                        if 'Protein' in reaction.id:
                            ignore_water = True
                        else:
                            ignore_water = False
                        res = get_sum_of_reaction(reaction, abs(biomass_reaction.metabolites[reactant]), ignore_water,
                                                  element_counter)
                        print(res)
                        c += res[0]
                        element_counter = res[1]
    return round(c, 3), element_counter


def get_light_kinetics(biomass, Eo=None, Lr=None, Ke=None):
    if not biomass: biomass = 0
    if not Eo:
        Eo = 1200
    if not Lr:
        Lr = 1
    if not Ke:
        Keo = 11.5
        wchl = 6 * 10e-3
        Ke = Keo * wchl * biomass
    E = Eo / (Lr * Ke) * (1 - np.exp(-Lr * Ke))
    return (0, E)


def get_micmen_kinetics(S, parameters):
    res = -parameters["Vmax"] * S / (parameters["Km"] + S)
    return res, 1000


def get_caro_kinetics(biomass, s=None, parameters=None):
    if parameters is None:
        parameters = {'wn': 0.03, 'n': 2}
    E = get_light_kinetics(biomass, s)[1]
    vcar_gen = car_gen(E, parameters['n'])
    a0 = 6.5e-2
    a1 = 7e-3 / 3600
    x = a1 * E + a0
    phi_val = phi(x - parameters['wn'])
    vcar = vcar_gen * phi_val
    return vcar, 1000


def phi(x):
    ns = 40
    return 1 / (1 + np.exp(-ns * x))


def car_gen(E, n):
    vmax = 8e-3 * 24
    E2 = E ** n
    Exa2 = (420 / 1000 * 24) ** n
    return vmax * E2 / (E2 + Exa2)


def nitrogen_quota():
    # TODO implement
    pass


def convert_mg_gDW_to_mmol_gDW(mg, MW):
    return mg / MW * 1000


def convert_mmol_mol_to_g_molMM(mmol, MW):
    return mmol * MW


def normalize(mg_gDW):
    return {key: mg_gDW[key] / sum(mg_gDW.values()) for key in mg_gDW.keys()}


def convert_mg_molMM_to_mmolM_gMM(mmol_molMM: dict, total):
    return {key: mmol_molMM[key] / total * 1000 for key in mmol_molMM.keys()}


def convert_mg_gMM_to_mmol_gMM(mg_gMM, MW):
    return mg_gMM / MW


# Disable
def block_print():
    sys.stdout = open(os.devnull, 'w')


def enable_print():
    sys.stdout = sys.__stdout__


def flux_change(fluxes_control: dict, fluxes_condition: dict, threshold: float = 0.1) -> dict:
    """
    Function to calculate the flux change between two conditions
    :param fluxes_control:
    :param fluxes_condition:
    :return:
    """
    flux_change = {}
    for key, value in fluxes_control.items():
        if abs(value + fluxes_condition[key]) != 0:
            flux_change[key] = (value - fluxes_condition[key]) / abs(value + fluxes_condition[key])
    fx = pd.DataFrame.from_dict(data=flux_change, orient='index', columns=['Flux_change'])
    as_dict = fx.loc[(fx.Flux_change >= threshold) | (fx.Flux_change <= -threshold)].to_dict(orient='index')
    for key, value in as_dict.items():
        as_dict[key] = round(value['Flux_change'], 3)
    return as_dict


def reaction_capacity(fva_solution: pd.DataFrame):
    """
    Function to calculate the reaction capacity
    :param fva_solution:
    :return:
    """
    as_dict = {}
    for index, row in fva_solution.iterrows():
        as_dict[index] = row['maximum'] - row['minimum']
    return as_dict


def differential_reaction_capacity(rc_1, rc_2):
    transformed_fold_changes = {}
    for key, value in rc_1.items():
        transformed_fold_changes[key] = math.log2(rc_2[key]) - math.log2(value)
    return transformed_fold_changes


def mp(function: callable, iterable, processes: int = 1, **kwargs):
    return Parallel(n_jobs=processes, **kwargs)(delayed(function)(iterable[0]))


def get_parameter_range(start, end, number_of_steps):
    linsp = np.linspace(round(start * 100, 0), round(end * 100, 0), int(round(number_of_steps * 100, 0)))
    return [e / 100 for e in linsp]
