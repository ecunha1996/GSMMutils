import subprocess

import numpy as np
from cobra import Metabolite
from scipy.stats import linregress
from experimental.BiomassComponent import BiomassComponent
from model.COBRAmodel import get_biomass_mass, MyModel


def get_productivity(data, time=48):
    previous_dw = 0
    productivity_g_l = [0]
    for i, row in data.iterrows():
        if int(i) == 0:
            previous_dw = row["DW"]
        if int(i) > 0:
            productivity_g_l.append((row["DW"] - previous_dw)/time)
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
            avv = (row["DW"] + previous_dw)/2
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


def get_element_in_biomass(model, element):
    res = get_biomass_mass(model)
    percentage_map = {}
    for key in res[1]:
        percentage_map[key] = round(res[1][key] * Metabolite(formula=key).formula_weight / 1000, 3)
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
    return round((b - a) / (exponential_phase[1] - exponential_phase[0]),3)

def get_maximum_productivity(data, exponential_phase):
    a = data.loc[str(exponential_phase[0])]["DW"]
    b = data.loc[str(exponential_phase[1])]["DW"]
    return round((b - a) / (exponential_phase[1] - exponential_phase[0]),3)


def run(cmd):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode("utf-8"))
    if stderr:
        print("\n\nError: ", stderr.decode("utf-8"))


def get_precursors(macromolecule, temp_precursor, model, biomass_dict):
    to_ignore = ["C00001__cytop"]
    precedent_reaction = [reaction for reaction in temp_precursor.reactions if "Biomass" not in reaction.id and
                          reaction.id.startswith("e_") and reaction.metabolites[temp_precursor] > 0]
    if len(precedent_reaction) > 1:
        raise Exception("More than one precedent reaction")
    if len(precedent_reaction) == 0:
        return temp_precursor.id
    else:
        if macromolecule.id not in to_ignore:
            for precursor in precedent_reaction[0].reactants:
                precursor = get_precursors(macromolecule, precursor, model, biomass_dict)
                if macromolecule.id in biomass_dict.keys() and type(precursor) == str:
                    precursor = BiomassComponent(precursor, precedent_reaction[0].metabolites[model.metabolites.get_by_id(precursor)], temp_precursor.id)
                    if temp_precursor == macromolecule:
                        biomass_dict[macromolecule.id].append(precursor)
                    else:
                        if temp_precursor.id not in biomass_dict[macromolecule.id]:
                            if type(biomass_dict[macromolecule.id]) == list:
                                biomass_dict[macromolecule.id] = {temp_precursor.id: [precursor]}
                            else:
                                biomass_dict[macromolecule.id][temp_precursor.id] = [precursor]
                        else:
                            biomass_dict[macromolecule.id][temp_precursor.id].append(precursor)
        return biomass_dict

def infer_biomass_from_model(model: MyModel, biomass_reaction_name):
    biomass_dict = {}
    biomass_reaction = model.reactions.get_by_id(biomass_reaction_name)
    macromolecules = biomass_reaction.reactants
    macromolecules = set(macromolecules) - {model.metabolites.get_by_id("C00001__cytop"), model.metabolites.get_by_id("C00002__cytop")}
    for macromolecule in macromolecules:
        biomass_dict[macromolecule.id] = []
        biomass_dict.update(get_precursors(macromolecule, macromolecule, model, biomass_dict))
    return biomass_dict