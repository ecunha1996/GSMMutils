import subprocess

import numpy as np
from cobra import Metabolite
from scipy.stats import linregress
from experimental.BiomassComponent import BiomassComponent


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


def get_precursors(macromolecule, temp_precursor, model):
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
                precursor = get_precursors(macromolecule, precursor, model)
                if type(precursor) == str:
                    if temp_precursor.id  in model.biomass_components:
                        parent_component = model.biomass_components[temp_precursor.id]
                    else:
                        parent_component = BiomassComponent(temp_precursor, precedent_reaction[0].metabolites[model.metabolites.get_by_id(precursor)], model.biomass_components[macromolecule.id])
                        model.biomass_components[temp_precursor.id] = parent_component
                    BiomassComponent(model.metabolites.get_by_id(precursor), precedent_reaction[0].metabolites[model.metabolites.get_by_id(precursor)], parent_component)


def update_st(stoichiometries, new_value):
    remaining = 0
    old_sum = 1 + sum(value for key, value in stoichiometries.items() if key in new_value.keys())
    new_sum = 1 + sum(new_value.values())
    for key, value in new_value.items():
        stoichiometries[key] = value
    for key, value in stoichiometries.items():
        if key not in new_value.keys():
            stoichiometries[key] = value * new_sum / old_sum
    print(sum(stoichiometries.values()))
    return stoichiometries


# def update_st(stoichiometries, new_value):
#     remaining = 0
#     new_value = {key: abs(value) for key, value in new_value.items()}
#     for key in new_value.keys():
#         remaining = remaining + abs(stoichiometries[key]) - new_value[key]
#         stoichiometries[key] = -new_value[key]
#     counter = 0
#     for key, value in stoichiometries.items():
#         if not key.startswith("C0"):
#             if key not in new_value.keys():
#                 counter += 1
#     for key, value in stoichiometries.items():
#         if not key.startswith("C0"):
#             if key not in new_value.keys():
#                 stoichiometries[key] += remaining / counter
#     print(sum(stoichiometries.values()))
#     return stoichiometries

def get_biomass_mass(model):
    def get_sum_of_reaction(reaction, stoichiometry, ignore_water,elementar_counter):
        counter = 0
        for reactant in reaction.reactants:
            copy = reactant.copy()
            copy.elements = {key: value for key, value in copy.elements.items() if key != "T"}
            counter += abs(reaction.metabolites[reactant]) * copy.formula_weight * stoichiometry
            for key in elementar_counter.keys():
                if key in reactant.elements:
                    elementar_counter[key] += abs(reaction.metabolites[reactant]) * reactant.elements[key] * stoichiometry
        for product in reaction.products:
            if product.id != "C00001__cytop":
                copy = product.copy()
                copy.elements = {key: value for key, value in copy.elements.items() if key != "T"}
                counter -= reaction.metabolites[product] * copy.formula_weight * stoichiometry
                for key in elementar_counter.keys():
                    if key in product.elements:
                        elementar_counter[key] -= abs(reaction.metabolites[product]) * product.elements[key] * stoichiometry
            else:
                if not ignore_water:
                    copy = product.copy()
                    copy.elements = {key: value for key, value in copy.elements.items() if key != "T"}
                    counter -= reaction.metabolites[product] * copy.formula_weight * stoichiometry
                    for key in elementar_counter.keys():
                        if key in product.elements:
                            elementar_counter[key] -= abs(reaction.metabolites[product]) * product.elements[key] * stoichiometry
        return counter / 1000, elementar_counter

    counter = 0
    elementar_counter = {"C":0,"N":0, "O":0, "P":0, "S":0, "H":0}
    lipid_stoichiometry = abs(model.reactions.e_Biomass__cytop.metabolites[model.metabolites.e_Lipid__cytop])
    for reactant in model.reactions.e_TAG__lip.reactants:
        counter += abs(model.reactions.e_TAG__lip.metabolites[reactant]) * reactant.formula_weight * 0.13 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_TAG__lip.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_DAG__er.reactants:
        counter += abs(model.reactions.e_DAG__er.metabolites[reactant]) * reactant.formula_weight * 0.1418 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_DAG__er.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_DGTS__er.reactants:
        counter += abs(model.reactions.e_DGTS__er.metabolites[reactant]) * reactant.formula_weight * 0.0863 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_DGTS__er.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_PE__er.reactants:
        counter += abs(model.reactions.e_PE__er.metabolites[reactant]) * reactant.formula_weight * 0.0432 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_PE__er.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_PC__er.reactants:
        counter += abs(model.reactions.e_PC__er.metabolites[reactant]) * reactant.formula_weight * 0.0386 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_PC__er.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_PI__er.reactants:
        counter += abs(model.reactions.e_PI__er.metabolites[reactant]) * reactant.formula_weight * 0.0329 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_PI__er.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_PG__chlo.reactants:
        counter += abs(model.reactions.e_PG__chlo.metabolites[reactant]) * reactant.formula_weight * 0.1808 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_PG__chlo.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_DGDG__chlo.reactants:
        counter += abs(model.reactions.e_DGDG__chlo.metabolites[reactant]) * reactant.formula_weight * 0.191 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_DGDG__chlo.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_SQDG__chlo.reactants:
        counter += abs(model.reactions.e_SQDG__chlo.metabolites[reactant]) * reactant.formula_weight * 0.114 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_SQDG__chlo.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_MGDG__chlo.reactants:
        counter += abs(model.reactions.e_MGDG__chlo.metabolites[reactant]) * reactant.formula_weight * 0.2373 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_MGDG__chlo.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_CL__mito.reactants:
        counter += abs(model.reactions.e_CL__mito.metabolites[reactant]) * reactant.formula_weight * 0.0264 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_CL__mito.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    for reactant in model.reactions.e_FFA__cytop.reactants:
        counter += abs(model.reactions.e_FFA__cytop.metabolites[reactant]) * reactant.formula_weight * 0.1153 * lipid_stoichiometry
        for key in elementar_counter.keys():
            if key in reactant.elements:
                elementar_counter[key] += abs(model.reactions.e_FFA__cytop.metabolites[reactant]) * reactant.elements[key] * 0.13 * lipid_stoichiometry
    c = 0
    c += counter / 1000
    to_ignore = ["C00002__cytop", "C00001__cytop", "e_Lipid__cytop"]
    for reactant in model.reactions.e_Biomass__cytop.reactants:
        if reactant.id not in to_ignore:
            for reaction in reactant.reactions:
                if reaction.id != "e_Biomass__cytop":
                    if 'Protein' in reaction.id:
                        ignore_water = True
                    else:
                        ignore_water = False
                    res = get_sum_of_reaction(reaction, abs(model.reactions.e_Biomass__cytop.metabolites[reactant]), ignore_water, elementar_counter)
                    c += res[0]
                    elementar_counter = res[1]
    return round(c,3), elementar_counter

