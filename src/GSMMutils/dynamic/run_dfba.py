import json
import numbers
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import shutil
import sys

from matplotlib import pyplot as plt

sys.path.insert(0, "/home/src/")
import pandas as pd
from parallelbar import progress_imap, progress_map
from timeout_decorator import timeout

from dfba import DfbaModel, ExchangeFlux, KineticVariable
import seaborn as sns
from functools import partial
from os.path import join
from GSMMutils.dynamic.initial_conditions import get_initial_conditions
from GSMMutils.experimental.ExpMatrix import *
from GSMMutils.model.COBRAmodel import MyModel
from GSMMutils.dynamic.rhs import get_bounds
from GSMMutils.dynamic.drhs import get_dynamic_expressions
from GSMMutils.dynamic.exchange_fluxes import get_exchange_fluxes
import sympy as sp
from sympy import Max, Min
from sympy.functions import Abs
from GSMMutils.graphics.plot import plot_concentrations
from tqdm import tqdm

from GSMMutils import DATA_PATH

os.chdir("/home/")

import random

matrix = ExpMatrix(f"{DATA_PATH}/experimental/Matriz- DCCR Dunaliella salina_dfba.xlsx")
matrix.conditions = "Resume"


def select_random_conditions(conditions, num_to_select, mandatory_strings):
    """
    Selects a specified number of strings from a list, making some of them mandatory.

    Args:
        conditions (list): List of strings to select from.
        num_to_select (int): Number of strings to randomly select.
        mandatory_strings (list): List of strings that must be included in the selected strings.

    Returns:
        list: List of selected strings.
    """
    # Check if the number of mandatory strings exceeds the number to select
    if len(mandatory_strings) > num_to_select:
        raise ValueError("Number of mandatory strings cannot exceed number of strings to select.")
    to_ignore = ["fachet", "Xi", "Yimei"]
    conditions = [e for e in conditions if not any(e.startswith(ignore) for ignore in to_ignore)]
    # Shuffle the input list of strings randomly
    random.shuffle(conditions)

    # Select the mandatory strings
    selected = mandatory_strings[:]

    # Select the remaining strings randomly
    remaining = num_to_select - len(mandatory_strings)
    selected.extend(random.sample(conditions, remaining))
    return selected


def read_model() -> MyModel:
    """
    Reads the model from the xml file and sets the objective function to the biomass reaction, demands, and minimizes the nitrate and phosphorus uptake.
    Returns
        MyModel: The model with the objective function set to the biomass reaction and the nitrate and phosphorus uptake minimized.
    -------

    """
    fba_model = MyModel(join(DATA_PATH, "models/model_dfba.xml"), "e_ActiveBiomass__cytop")
    fba_model.exchanges.EX_C00011__dra.bounds = (-10000, 10000)
    fba_model.solver = "glpk"
    [setattr(x, 'objective_coefficient', 0) for x in fba_model.reactions if x.objective_coefficient != 0]

    objectives = {"e_ActiveBiomass__cytop": 1, "DM_C00369__chlo": 1, "DM_C05306__chlo": 1, "DM_C05307__chlo": 1, "DM_C08601__chlo": 1,
                  "DM_C02094__chlo": 1, "DM_C00116__cytop": 1, "DM_C00422__lip": 1, "EX_C00244__dra": -1, "EX_C00009__dra": -1
                  }
    for reaction_id, value in objectives.items():
        fba_model.reactions.get_by_id(reaction_id).objective_coefficient = value
    return fba_model


fba_model = read_model()


def get_kinetic_variables() -> dict:
    """
    Creates the kinetic variables used in the model.
    Returns
    -------

    """
    light = KineticVariable("Light")
    X = KineticVariable("Biomass")
    F = KineticVariable("ActiveBiomass")
    nitrate = KineticVariable("Nitrate")
    n_quota = KineticVariable("Nitrogen_quota")
    phosphorus = KineticVariable("Phosphate")
    p_quota = KineticVariable("Phosphate_quota")
    starch = KineticVariable("Starch")
    starch_concentration = KineticVariable("Starch_concentration")
    tag = KineticVariable("TAG")
    glycerol = KineticVariable("Glycerol")
    carotene = KineticVariable("Carotene")
    lutein = KineticVariable("Lutein")
    chlorophyll = KineticVariable("Chlorophyll")
    return {"X": X, "F": F, "nitrate": nitrate, "phosphorus": phosphorus, "starch": starch, "starch_concentration": starch_concentration, "chlorophyll": chlorophyll, "carotene": carotene,
            "n_quota": n_quota, "p_quota": p_quota, "tag": tag, "glycerol": glycerol, "light": light, "lutein": lutein}


@timeout(20)
def create_dfba_model(condition, parameters, create_plots=False):
    """
    Creates the dfba model for a given condition and parameters.
    Parameters
    ----------
    condition (str): Condition to create the model for.
    parameters (dict): Dictionary with the parameters to use in the model.
    create_plots (bool): Whether to create plots or not.

    Returns
    -------

    """
    if 'Time (d)' not in matrix.matrix[condition].columns:
        matrix.matrix[condition]['Time (d)'] = matrix.matrix[condition].index
    dfba_model = DfbaModel(fba_model)
    dfba_model.solver_data.set_display("none")
    dfba_model.solver_data.set_algorithm("direct")
    #     dfba_model.solver_data.set_rel_tolerance(1e-3)
    #     dfba_model.solver_data.set_ode_method("ADAMS")

    kinetic_vars = get_kinetic_variables()
    parameters.update(kinetic_vars)
    dfba_model.add_kinetic_variables(list(kinetic_vars.values()))

    exchange_fluxes = {}
    for key, value in get_exchange_fluxes().items():
        exchange_fluxes[key] = ExchangeFlux(value)

    dfba_model.add_exchange_fluxes(list(exchange_fluxes.values()))
    parameters.update(exchange_fluxes)
    parameters['starch_production'] = parameters['v_S'] * 48660.195 / 1000
    parameters['chl_production'] = parameters['v_chla'] * 893.49 / 1000 + parameters['v_chlb'] * 907.49 / 1000
    parameters['caro_production'] = parameters['v_C'] * 536.87 / 1000
    parameters['lutein_production'] = parameters['v_lutein'] * 568.87 / 1000
    parameters['glycerol_production'] = parameters['v_glycerol'] * 92.09 / 1000
    parameters['tag_production'] = parameters['v_tag'] * 904.78 / 1000
    parameters['total_growth_rate'] = parameters['mu'] + parameters['starch_production'] + parameters['chl_production'] + parameters['caro_production'] + parameters['glycerol_production'] + parameters['tag_production'] + parameters['lutein_production']

    for key, value in get_dynamic_expressions(parameters).items():
        dfba_model.add_rhs_expression(key, value)


    ##### Experiment-dependent ######
    parameters["Eo"] = matrix.conditions["Light (umol/m^2.s)"].loc[condition]
    parameters["nacl"] = matrix.conditions["Salinity g/L"].loc[condition]
    parameters["Lr"] = matrix.conditions["Lr"].loc[condition]
    parameters["aeration"] = matrix.conditions["Aeration rate"].loc[condition]

    light_sources = matrix.conditions["Light sources"].loc[condition].split(",")
    for light_source in light_sources:
        fba_model.reactions.get_by_id(light_source).bounds = (0, 10000)

    ##### General ######
    parameters["q"] = parameters["n_quota"] / parameters["wNmax"]
    parameters["n"] = 1 - (parameters["q"] / (parameters["q"] + parameters["K_nitrogen_quota"]))
    x_storage = parameters["starch"] + parameters["carotene"] + parameters["glycerol"] + parameters["tag"]
    cell_size_increase = 1 / (1 - x_storage)
    parameters["z"] = (cell_size_increase - 1) / (parameters["t_max"] - 1)
    parameters["nitrogen_mass_quota"] = parameters["n_quota"] * 14.01 / 1000
    parameters["phosphate_mass_quota"] = parameters["p_quota"] * 30.97 / 1000

    ##### Light ######
    parameters["Ex"] = get_bounds("light", parameters)
    dfba_model.add_exchange_flux_lb("EX_C00205__dra", sp.Max(sp.N(parameters["Ex"]), 0), parameters["light"])

    ##### NO3 ######
    dfba_model.add_exchange_flux_lb("EX_C00244__dra", get_bounds("nitrate", parameters), parameters["nitrate"])  # 4.07
    #     nitrate_quota = sp.Max(0, 1 - (4.8697 * F / X) / n_quota)
    dfba_model.add_exchange_flux_lb("DM_C00244__cytop", sp.Max(0, parameters["v_nitrate_max"] * (1 - parameters['wNmin'] / parameters["n_quota"])), parameters["nitrate"])

    ##### HPO4 ######
    dfba_model.add_exchange_flux_lb("EX_C00009__dra", sp.Max(0, parameters["VPmax"] * parameters["phosphorus"] / (parameters['KPm'] + parameters["phosphorus"])), parameters["phosphorus"])
    #     polyP_quota = sp.Max(0, 1 - (0.295 * F / X) / p_quota)
    dfba_model.add_exchange_flux_lb("DM_C00404__vacu", sp.Max(0, parameters["v_polyphosphate_max"] * (1 - parameters["wPmin"] / parameters["p_quota"])), parameters["phosphorus"])

    ##### Starch ######
    dfba_model.add_exchange_flux_lb("DM_C00369__chlo", get_bounds("starch_consumption", parameters), parameters["starch"])
    dfba_model.add_exchange_flux_ub("DM_C00369__chlo", get_bounds("starch_production", parameters), parameters["starch"])  #

    ##### Carotene ######
    dfba_model.add_exchange_flux_ub("DM_C02094__chlo", get_bounds("carotene", parameters), parameters["carotene"])

    ##### Lutein ######
    dfba_model.add_exchange_flux_ub("DM_C08601__chlo", get_bounds("lutein", parameters), parameters["lutein"])

    ##### Chlorophyll ######
    sum_chl = get_bounds("chlorophyll", parameters)
    dfba_model.add_exchange_flux_lb("DM_C05306__chlo", Abs(Min(sum_chl, 0)) * 1.73 / 2.73, parameters["chlorophyll"])
    dfba_model.add_exchange_flux_lb("DM_C05307__chlo", Abs(Min(sum_chl, 0)) / 2.73, parameters["chlorophyll"])
    dfba_model.add_exchange_flux_ub("DM_C05306__chlo", Max(sum_chl, 0) * 1.73 / 2.73, parameters["chlorophyll"])
    dfba_model.add_exchange_flux_ub("DM_C05307__chlo", Max(sum_chl, 0) / 2.73, parameters["chlorophyll"])

    ##### Glycerol ######
    # wgly_max = 0.17  # https://doi.org/10.1016/j.biortech.2008.02.042
    dfba_model.add_exchange_flux_ub("DM_C00116__cytop", get_bounds("glycerol", parameters), parameters["glycerol"])

    ##### TAG ######
    dfba_model.add_exchange_flux_ub("DM_C00422__lip", get_bounds("tag", parameters), parameters["tag"])

    ##### CO2 ######
    dfba_model.add_exchange_flux_lb("EX_C00011__dra", get_bounds("co2", parameters))  # vco2max * (1 - z)
    dfba_model.add_initial_conditions(
        get_initial_conditions(matrix, condition)
    )
    max_time = max(matrix.matrix[condition]['Time (d)'].astype(float).tolist()) + 1
    time_step = 1 / 48
    concentrations, trajectories = dfba_model.simulate(
        0.0, max_time, time_step, ["e_ActiveBiomass__cytop", 'EX_C00009__dra', 'DM_C00369__chlo', "DM_C02094__chlo", "EX_C00244__dra", 'DM_C05306__chlo', "DM_C05307__chlo", "EX_C00011__dra",
                                   "DM_C00116__cytop", "DM_C00404__vacu", "DM_C00244__cytop", "EX_C00205__dra", "DM_C08601__chlo"]
    )

    active_biomass_fraction = concentrations['ActiveBiomass'] / concentrations['Biomass']
    concentrations['Protein'] = abs(fba_model.reactions.e_ActiveBiomass__cytop.metabolites[fba_model.metabolites.e_Protein__cytop]) * active_biomass_fraction
    carbs = abs(fba_model.reactions.e_ActiveBiomass__cytop.metabolites[fba_model.metabolites.e_Carbohydrate__cytop]) * active_biomass_fraction
    concentrations['Carbohydrate'] = carbs + concentrations['Starch']
    polar_lipids = abs(fba_model.reactions.e_ActiveBiomass__cytop.metabolites[fba_model.metabolites.e_Lipid__cytop]) * active_biomass_fraction
    concentrations['Lipid'] = polar_lipids + concentrations['TAG']

    indexes = matrix.matrix[condition].index.astype(float)
    indexes = [e for e in indexes]
    experimental = [(indexes, matrix.matrix[condition]["DW"].tolist())]

    concentrations['Carotene_concentration'] = concentrations['Carotene'] * concentrations['Biomass']
    concentrations['Chlorophyll_concentration'] = concentrations['Chlorophyll'] * concentrations['Biomass']
    concentrations['Lutein_concentration'] = concentrations['Lutein'] * concentrations['Biomass']

    if create_plots:

        if 'Caro' in matrix.matrix[condition].columns:
            if condition.startswith("fachet") or condition.startswith("Xi") or condition.startswith("Yimei"):
                molecules = ['Caro']
                experimental_caro = [[matrix.matrix[condition][molecule].dropna().index.astype(float).tolist(), matrix.matrix[condition][molecule].dropna().tolist()] for molecule in molecules]
                plot_concentrations(concentrations, y=['Carotene'], experimental=experimental_caro, filename=f"{DATA_PATH}/dfba/pigments/carotene_{condition}.png", y_label="Macromolecule (g/gDW)", experimental_label=molecules)

        if 'Chl' in matrix.matrix[condition].columns:
            if condition.startswith("fachet") or condition.startswith("Xi") or condition.startswith("Yimei"):
                molecules = ['Chl']
                experimental_caro = [[matrix.matrix[condition][molecule].dropna().index.astype(float).tolist(), matrix.matrix[condition][molecule].dropna().tolist()] for molecule in molecules]
                plot_concentrations(concentrations, y=['Chlorophyll'], experimental=experimental_caro, filename=f"{DATA_PATH}/dfba/pigments/chlorophyll_{condition}.png", y_label="Macromolecule (g/gDW)", experimental_label=molecules)

        if 'Caro_concentration' in matrix.matrix[condition].columns:
            if condition.startswith("fachet") or condition.startswith("Xi") or condition.startswith("Yimei"):
                molecules = ['Caro_concentration']
                experimental_caro = [[matrix.matrix[condition][molecule].dropna().index.astype(float).tolist(), matrix.matrix[condition][molecule].dropna().tolist()] for molecule in molecules]
                plot_concentrations(concentrations, y=['Carotene_concentration'], experimental=experimental_caro, filename=f"{DATA_PATH}/dfba/pigments/carotene_conc_{condition}.png", y_label="Macromolecule (g/L)", experimental_label=molecules)

        if 'Chlorophyll_concentration' in matrix.matrix[condition].columns:
            if condition.startswith("fachet") or condition.startswith("Xi") or condition.startswith("Yimei"):
                molecules = ['Chlorophyll_concentration']
                experimental_caro = [[matrix.matrix[condition][molecule].dropna().index.astype(float).tolist(), matrix.matrix[condition][molecule].dropna().tolist()] for molecule in molecules]
                plot_concentrations(concentrations, y=['Chlorophyll_concentration'], experimental=experimental_caro, filename=f"{DATA_PATH}/dfba/pigments/chlorophyll_conc_{condition}.png", y_label="Macromolecule (g/L)", experimental_label=molecules)

        if 'NO3' in matrix.matrix[condition].columns:
            molecules = ['NO3']
            experimental_no3 = [[matrix.matrix[condition][molecule].dropna().index.astype(float).tolist(), matrix.matrix[condition][molecule].dropna().tolist()] for molecule in molecules]
            plot_concentrations(concentrations, y=['Nitrate'], experimental=experimental_no3, filename=f"{DATA_PATH}/dfba/nitrate_{condition}.png", y_label="Macromolecule (g/gDW)", experimental_label=molecules)

        fig1 = plot_concentrations(concentrations, y=["Biomass", "ActiveBiomass"], experimental=experimental, filename=f"{DATA_PATH}/dfba/biomass_concentrations/biomass_concentrations_{condition}.png", y_label="Biomass (g/L)")

        fig2 = plot_concentrations(concentrations, y=["Phosphate"], secondary_axis=["Nitrate"], filename=f"{DATA_PATH}/dfba/concentrations/external_concentrations_{condition}.png", y_label="Phosphate (mmol/L)", secondary_y_label="Nitrate (mmol/L)")

        fig3 = plot_concentrations(concentrations, y=["Glycerol", "Starch", "TAG"], filename=f"{DATA_PATH}/dfba/quotas/intracellular_quotas_{condition}.png", y_label="Quota (g/gDW)", secondary_axis=["Carotene", "Chlorophyll"], secondary_y_label="Quota (g/gDW)")

        molecules = ["Protein", "Lipid", "Carbohydrate"]
        if all(molecule in matrix.matrix[condition].columns for molecule in molecules):
            experimental = [[matrix.matrix[condition][molecule].dropna().index.astype(float).tolist(), matrix.matrix[condition][molecule].dropna().tolist()] for molecule in molecules]
            fig4 = plot_concentrations(concentrations, y=molecules, experimental=experimental, filename=f"{DATA_PATH}/dfba/macros/macromolecules_{condition}.png", y_label="Macromolecule (g/gDW)", experimental_label=molecules)

        concentrations.to_csv(f"{DATA_PATH}/dfba/concentrations/concentrations_{condition}.csv", index=False)
        trajectories.to_csv(f"{DATA_PATH}/dfba/trajectories/trajectories_{condition}.csv", index=False)

    return concentrations, trajectories


def get_closest(list_a, list_b):
    """
    Gets the closest values from list_b to list_a.
    Parameters
    ----------
    list_a
    list_b

    Returns
    -------

    """
    closest_values = []
    for num_A in list_a:
        index = np.abs(np.array(list_b) - num_A).argmin()
        closest_values.append(list_b[index])
    return closest_values


def fitness(conditions_names, parameters_under_optimization=None, initial_parameters=None):
    """
    Calculates the fitness of a set of parameters.
    Parameters
    ----------
    conditions_names
    parameters_under_optimization
    initial_parameters

    Returns
    -------

    """
    if not parameters_under_optimization:
        parameters = {}
        parameters_names = ('ro1', 'ro0', 'a0', 'a1', 'a2', 'a3', 'a4', 'ExA', 'l', 'smoothing_factor', 'wPopt', 'wPmin', 'wNmax', 'wNmin', 'c0', 't_max', 'K_nitrogen_quota', 'VPmax', 'KPm', 'VNmax', 'KNm', 'wgly_max',
                            'maximum_starch_production', 'maximum_tag_production', 'v_nitrate_max', 'v_polyphosphate_max', "v_car_max", "ymax", "Esat", "KEchl", "vco2max", "Kstl", "hill_coeff_starch", "light_conversion_factor",
                            "v_lut_max", "a0_lut", "a1_lut", "a3_lut", "a4_lut", "smoothing_factor_lut", "nacl_lipid")

        for i, parameter_name in enumerate(parameters_names):
            parameters[parameter_name] = 2 ** initial_parameters[i]
    else:
        parameters = json.load(open(f"{DATA_PATH}/dfba/inputs/initial_parameters.json", "r"))
        for index in range(len(initial_parameters)):
            parameters[parameters_under_optimization[index]] = 2 ** initial_parameters[index]
    # with open(f"{DATA_PATH}/dfba/inputs/parameters.json", "w") as f:
    #     json.dump(parameters, f, indent=4)
    # print(f"Parameters: {parameters}")
    try:
        total_error = sum(Parallel(n_jobs=len(conditions_names), timeout=30, backend="multiprocessing")(delayed(evaluate_trial)(parameters, condition=condition) for condition in conditions_names))
        # total_error = sum(progress_imap(partial(evaluate_trial, parameters, True), conditions_names,
        #                                 process_timeout = 30))
        # total_error = 0
        # pbar = tqdm(total=len(conditions_names), desc="Evaluating trials")
        # for condition in conditions_names:
        #     total_error += evaluate_trial(parameters, True, condition)
        #     pbar.set_description(f"Evaluating {condition}; error is {total_error}")
        #     pbar.update(1)
    except Exception as e:
        print(e)
        total_error = 1e3
    print(f"Total error from set of parameters: {total_error}")
    with open(f"{DATA_PATH}/dfba/temp_error.log", "w") as file:
        file.write(f"{total_error}\n")
    return round(total_error, 2)


def evaluate_trial(parameters, create_plots=False, condition=None):
    """
    Evaluates a trial.
    Parameters
    ----------
    parameters (dict): Dictionary with the parameters to use in the model.
    create_plots (bool): Whether to create plots or not.
    condition (str): Condition to create the model for.

    Returns (float): The total error of the trial.
    -------

    """
    print(f"Trial: {condition}")
    total_error = 0
    mat = matrix.matrix[condition]
    mat['Time (d)'] = [round(e, 2) for e in mat.index.astype(float)]
    try:
        concentrations, trajectories = create_dfba_model(condition, parameters, create_plots)
        # to_fit = {"Biomass": "DW", "Carotene": "Caro", "Chlorophyll": "Chl", "Starch": "Starch", "Nitrate": "NO3", 'Protein': 'Protein', 'Carbohydrate': 'Carbohydrate', 'Lipid': 'Lipid',
        #           "Chlorophyll_concentration": "Chlorophyll_concentration", "Carotene_concentration": "Caro_concentration"
        #           }  #
        to_fit = {"Biomass": "DW", 'Lipid': 'Lipid'}
        experimental_time = np.array(mat["Time (d)"])
        closest = get_closest(experimental_time, concentrations.time)
        at_time = concentrations.loc[concentrations.time.isin(closest)]
        at_time.reset_index(inplace=True, drop=True)
        mat.reset_index(inplace=True, drop=True)
        for simulation_name, experimental_name in to_fit.items():
            if experimental_name in mat.columns:
                experimental = mat[experimental_name]
                simulated = at_time[simulation_name]
                total_error += get_relative_error(experimental, simulated)
                print(f"Total error for {simulation_name}:\n{total_error}")
    except Exception as e:
        print(e)
        with open(f"{DATA_PATH}/dfba/temp_error.log", "a") as file:
            file.write(f"{e}\n")
        total_error = 1e3
    return round(total_error, 3)


def get_relative_error(experimental, simulated):
    """
    Calculates the relative error between two vectors.
    Parameters
    ----------
    experimental
    simulated

    Returns
    -------

    """
    relative_error = 0
    if simulated.shape[0] < experimental.shape[0] / 4:
        return 100
    try:
        experimental.dropna(inplace=True)
        intersection = [e for e in simulated.index if e in experimental.index]
        simulated = simulated.loc[intersection]
        experimental = experimental.loc[intersection]
        abs_error = np.abs(experimental - simulated)
        with np.errstate(divide='ignore', invalid='ignore'):
            relative_error = round(sum(np.where(experimental == 0, np.inf, abs_error / np.abs(experimental))), 4)
    #         print(f"Relative error:\n{relative_error}")
    except Exception as e:
        print(e)
        print(f"Experimental values:\n{experimental}")
        print(f"Simulated values:\n{simulated}")
    if not isinstance(relative_error, numbers.Number):
        print(f"Relative error:\n{relative_error}")
        print(f"Experimental values:\n{experimental}")
        print(f"Simulated values:\n{simulated}")
    if np.isnan(relative_error):
        print("NaN found!!")
        print(f"Experimental values:\n{experimental}")
        print(f"Simulated values:\n{simulated}")
        relative_error = 100
    return relative_error


def callbackF(pbar, x):
    """
    Callback function for the optimization.
    Parameters
    ----------
    pbar
    x

    Returns
    -------

    """
    with open(f"{DATA_PATH}/dfba/temp_error.log") as file:
        error = float(file.read())
    pbar.set_description(f"Current objective: {round(error, 3)}")
    pbar.update(1)


def parameter_optimization():
    """
    Runs the parameter optimization.
    Returns
    -------

    """
    from scipy.optimize import minimize

    initial_parameters = json.load(open(f"{DATA_PATH}/dfba/inputs/initial_parameters.json", "r"))

    # Define optimization method ('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', etc.)
    # Refer to scipy.optimize documentation for more options
    method = 'Nelder-Mead'

    # Define bounds for parameters (if any)
    bounds = [(0, 100), (0, 0.1), (0, 0.5), (0, 0.5), (0, 100), (0, 100), (0, 1), (0, 500), (0, 10), (0, 100), (0, 2), (0, 0.15), (0, 7), (0, 3), (0, 0.05), (1, 4), (0, 1), (0, 1), (0, 0.5), (0, 20), (0, 1), (0, 0.5), (0, 0.1), (0, 100), (0, 10),
              (0, 0.1), (0, 0.5), (0, 1), (100, 1500), (0, 150), (10, 150), (0, 2000), (0, 10), (0, 10), (0, 0.5), (0,1), (0,1), (0, 100), (0, 1), (0, 10), (0,0.1)]  # , (0, 10)
    conditions_names = set(matrix.matrix.keys()) - {"Resume"}
    validation = select_random_conditions(list(conditions_names), 5, ['fachet_HLND', "Yimei_HL"])
    conditions_names = tuple(conditions_names - set(validation) - {e for e in conditions_names if e.startswith("Xi") or e.startswith("Yimei")
                                                                   or e.startswith("fachet")})
    # conditions_names  = tuple(['Xi_cont_S', 'Xi_cont_F'])

    # pbar = tqdm(conditions_names, desc="Running initial conditions")
    # initial_error = sum(Parallel(n_jobs=len(conditions_names))(delayed(evaluate_trial)(initial_parameters, True, condition=condition) for condition in pbar))
    # pbar.set_description(f"Initial error: {initial_error}")
    initial_error = sum(progress_imap(partial(evaluate_trial, initial_parameters, True), conditions_names, n_cpu=len(conditions_names)))
    #
    #
    with open(f"{DATA_PATH}/dfba/validation.txt", 'w') as f:
        f.write("Validation conditions:\n")
        for e in validation:
            f.write(f"{e}\n")
        f.write(f"Initial error was: {initial_error}")
    shutil.make_archive(f'{DATA_PATH}/dfba', 'zip', f'{DATA_PATH}/dfba')
    max_iterations = 1000
    initial_parameters_log = [math.log2(e) for e in initial_parameters.values()]
    bounds_log = [(math.log2(e[0] + 1e-10), math.log2(e[1] + 1e-10)) for e in bounds]
    with tqdm(total=max_iterations, desc=f"Running optimization for {len(conditions_names)}") as pbar:
        result = minimize(partial(fitness, conditions_names, None), np.array(initial_parameters_log), method=method, bounds=bounds_log, callback=partial(callbackF, pbar),
                          options={"maxiter": max_iterations})
    # with tqdm(total=max_iterations, desc="Running optimization") as pbar:
    #     result = differential_evolution(partial(fitness, conditions_names, None), bounds_log, strategy='best1bin', maxiter=max_iterations, popsize=5, callback=partial(callbackF, pbar),
    #                                           workers=5, x0=initial_parameters_log)

    # Extract optimized parameters and fitness value
    optimal_params = [2 ** e for e in result.x]
    for index, param in enumerate(optimal_params):
        if round(param, 5) != 0:
            optimal_params[index] = round(param, 5)

    optimal_fitness = result.fun
    optimal_parameters = {list(initial_parameters.keys())[index]: param for index, param in enumerate(optimal_params)}
    # Print optimized parameters and fitness value
    print("Optimized Parameters:\n")

    with open(f"{DATA_PATH}/dfba/optimized_parameters.json", "w") as f:
        json.dump(optimal_parameters, f, indent=4)

    with open(f"{DATA_PATH}/dfba/optimized_parameters.txt", "w") as f:
        for index, param in enumerate(optimal_params):
            f.write(f"{list(initial_parameters.keys())[index]}\t{param}\n")
            print(f"{list(initial_parameters.keys())[index]}:\t{param}\t{optimal_params[index]}\n")
        f.write(f"Error after optimization:\t{str(optimal_fitness)}")
    print("Optimized Fitness Value: ", optimal_fitness)
    # final_error = sum(Parallel(n_jobs=30)(delayed(evaluate_trial)(fba_model, matrix, condition, optimal_parameters, create_plots=True) for condition in conditions_names))
    # print(f"Final error was: {final_error}")

    with tqdm(total=len(validation), desc="Running validation") as pbar:
        validation_error = sum(Parallel(n_jobs=30)(delayed(evaluate_trial)(optimal_parameters, create_plots=True, condition=condition) for condition in validation))
    pbar.set_description(f"Validation error: {validation_error}")
    # validation_error = sum(progress_imap(partial(evaluate_trial, fba_model, matrix, optimal_parameters, True), conditions_names, n_cpu=len(conditions_names)))
    print(f"Validation error was: {validation_error}")
    with open(f"{DATA_PATH}/dfba/validation.txt", 'a') as f:
        f.write(f"\nValidation error was: {validation_error}")


def optimize_simple_parameters():
    """
    Runs the parameter optimization.
    Returns
    -------

    """
    from scipy.optimize import minimize
    parameters_to_optimize = [
        'a0',
        'a1',
        'a2',
        'a3',
        'a4',
        'c0',
        'l',
        'light_conversion_factor',
        'ro0',
        'ro1',
        'smoothing_factor',
        'v_car_max',
    ]

    initial_parameters = json.load(open(f"{DATA_PATH}/dfba/inputs/initial_parameters.json", "r"))
    method = 'Nelder-Mead'
    bounds = [(0, 5), (0, 0.1), (0, 100), (0, 100), (0, 0.1), (0, 0.01), (0, 10), (0, 10), (0, 0.1), (0, 10), (0, 10), (0, 0.1)]
    conditions_names = set(matrix.matrix.keys()) - {"Resume"}
    conditions_names = set([condition for condition in conditions_names if not condition.startswith("fachet") and not condition.startswith("Xi") and not condition.startswith("Yimei")])
    validation = select_random_conditions(list(conditions_names), 5, ['fachet_ML'])
    conditions_names = tuple(conditions_names - set(validation))
    # conditions_names = tuple([condition for condition in conditions_names if not condition.startswith("fachet") and not condition.startswith("Xi") and not condition.startswith("Yimei")])
    initial_error = sum(progress_imap(partial(evaluate_trial, initial_parameters, True), conditions_names, n_cpu=len(conditions_names)))
    with open(f"{DATA_PATH}/dfba/validation.txt", 'w') as f:
        f.write("Validation conditions:\n")
        for e in validation:
            f.write(f"{e}\n")
        f.write(f"Initial error was: {initial_error}")
    shutil.make_archive(f'{DATA_PATH}/dfba', 'zip', f'{DATA_PATH}/dfba')
    max_iterations = 500
    initial_parameters = {k: v for k, v in initial_parameters.items() if k in parameters_to_optimize}
    initial_parameters_log = [math.log2(e) for e in initial_parameters.values()]
    bounds_log = [(math.log2(e[0] + 1e-10), math.log2(e[1] + 1e-10)) for e in bounds]
    with tqdm(total=max_iterations, desc=f"Running optimization for {len(conditions_names)}") as pbar:
        result = minimize(partial(fitness, conditions_names, parameters_to_optimize), np.array(initial_parameters_log), method=method, bounds=bounds_log, callback=partial(callbackF, pbar),
                          options={"maxiter": max_iterations})
    # with tqdm(total=max_iterations, desc="Running optimization") as pbar:
    #     result = differential_evolution(partial(fitness, conditions_names), bounds, strategy='best1bin', maxiter=max_iterations, popsize=5, callback=partial(callbackF, pbar),
    #                                           workers=2, x0=list(initial_parameters.values()))

    # Extract optimized parameters and fitness value
    optimal_params = [2 ** e for e in result.x]
    for index, param in enumerate(optimal_params):
        if round(param, 5) != 0:
            optimal_params[index] = round(param, 5)

    optimal_fitness = result.fun
    optimal_parameters = {list(initial_parameters.keys())[index]: param for index, param in enumerate(optimal_params)}
    # Print optimized parameters and fitness value
    print("Optimized Parameters:\n")

    with open(f"{DATA_PATH}/dfba/optimized_parameters.json", "w") as f:
        json.dump(optimal_parameters, f)

    with open(f"{DATA_PATH}/dfba/optimized_parameters.txt", "w") as f:
        for index, param in enumerate(optimal_params):
            f.write(f"{list(initial_parameters.keys())[index]}\t{param}\n")
            print(f"{list(initial_parameters.keys())[index]}:\t{param}\t{optimal_params[index]}\n")
        f.write(str(optimal_fitness))
    print("Optimized Fitness Value: ", optimal_fitness)
    # final_error = sum(Parallel(n_jobs=30)(delayed(evaluate_trial)(fba_model, matrix, condition, optimal_parameters, create_plots=True) for condition in conditions_names))
    # print(f"Final error was: {final_error}")

    with tqdm(total=len(validation), desc="Running validation") as pbar:
        validation_error = sum(Parallel(n_jobs=30)(delayed(evaluate_trial)(optimal_parameters, create_plots=True, condition=condition) for condition in validation))
    pbar.set_description(f"Validation error: {validation_error}")
    # validation_error = sum(progress_imap(partial(evaluate_trial, fba_model, matrix, optimal_parameters, True), conditions_names, n_cpu=len(conditions_names)))
    print(f"Validation error was: {validation_error}")
    with open(f"{DATA_PATH}/dfba/validation.txt", 'a') as f:
        f.write(f"\nValidation error was: {validation_error}")


def run_all_parallel():
    """
    Runs all the conditions in parallel.
    Returns
    -------

    """
    initial_parameters = json.load(open(f"{DATA_PATH}/dfba/inputs/initial_parameters.json", "r"))

    # model = read_model()
    # matrix = ExpMatrix(f"{DATA_PATH}/experimental/Matriz- DCCR Dunaliella salina_dfba.xlsx")
    # matrix.conditions = "Resume"
    conditions_names = sorted(tuple(set(matrix.matrix.keys()) - {"Resume"}))
    conditions_names = set([condition for condition in conditions_names if not condition.startswith("fachet") and not condition.startswith("Xi") and not condition.startswith("Yimei")])
    # total_error = sum(Parallel(n_jobs=30, timeout=500)(delayed(evaluate_trial)(model, matrix, initial_parameters, create_plots=True, condition=condition) for condition in conditions_names))
    total_error = sum(progress_map(partial(evaluate_trial, initial_parameters, True), conditions_names, n_cpu=len(conditions_names), process_timeout=30))
    with open(f"{DATA_PATH}/dfba/total_error.txt", "w") as f:
        f.write(str(total_error))
    st = {
        "1": "0.006112313",
        "2": "0.026350859",
        "3": "0.014449316",
        "4": "0.026350859",
        "5": "0.014449316",
        "6": "0.014789985",
        "7": "0.014449316",
        "8": "0.014789985",
        "9": "0.022663961",
        "10": "0.014789985",
        "11": "0.022663961",
        "12": "0.003563663",
        "13": "0.022663961",
        "14": "0.003563663",
        "15": "0.041744485",
        "16": "0.003563663",
        "17": "0.041744485",
        "18": "0.003208629",
        "20": "0.041744485",
        "22": "0.003208629",
        "23": "0.006232213",
        "24": "0.003208629",
        "RPC1": "0.00221518",
        "RPC2": "0.006232213",
        "RPC3": "0.00221518"
    }
    std_devs = {key: float(value) / 1000 for key, value in st.items()}
    error = list(std_devs.values())
    data_carotene, data_carotene_conc, data_chl, data_protein, data_lipid, data_carbohydrate, data_lutein, data_lipid_conc = {}, {}, {}, {}, {}, {}, {}, {}
    experimental_data_carotene, experimental_data_carotene_concentration, experimental_data_chl, experimental_data_protein, \
        experimental_data_lipid, experimental_data_carbohydrate, experimental_data_lutein, experimental_data_lipid_concentration= {}, {},{}, {}, {}, {}, {}, {}
    for condition in matrix.conditions.index:
        if 'Caro' in matrix.matrix[condition].columns and not condition.startswith("fachet") and not condition.startswith("Xi") and not condition.startswith("Yimei"):
            temp = pd.read_csv(f"{DATA_PATH}/dfba/concentrations/concentrations_{condition}.csv")
            data_carotene[condition] = temp.iloc[-1]['Carotene']
            data_carotene_conc[condition] = temp.iloc[-1]['Carotene'] * temp.iloc[-1]['Biomass']
            data_lutein[condition] = temp.iloc[-1]['Lutein']
            data_chl[condition] = temp.iloc[-1]['Chlorophyll'] * temp.iloc[-1]['Biomass']
            data_protein[condition] = temp.iloc[-1]['Protein']
            data_lipid[condition] = temp.iloc[-1]['Lipid']
            data_lipid_conc[condition] = temp.iloc[-1]['Lipid'] * temp.iloc[-1]['Biomass']
            data_carbohydrate[condition] = temp.iloc[-1]['Carbohydrate']
            experimental_data_carotene[condition] = matrix.matrix[condition]['Caro'].dropna().tolist()[-1]
            experimental_data_lutein[condition] = matrix.matrix[condition]['Lutein'].dropna().tolist()[-1]
            experimental_data_carotene_concentration[condition] = matrix.matrix[condition]['Caro'].dropna().tolist()[-1] * matrix.matrix[condition]['DW'].dropna().tolist()[-1]
            experimental_data_chl[condition] = matrix.matrix[condition]['Chl'].dropna().tolist()[-1] * matrix.matrix[condition]['DW'].dropna().tolist()[-1]
            molecules = ["Protein", "Lipid", "Carbohydrate"]
            if all(molecule in matrix.matrix[condition].columns for molecule in molecules):
                experimental_data_protein[condition] =  matrix.matrix[condition]['Protein'].dropna().tolist()[-1]
                experimental_data_lipid[condition] = matrix.matrix[condition]['Lipid'].dropna().tolist()[-1]
                experimental_data_lipid_concentration[condition] = matrix.matrix[condition]['Lipid'].dropna().tolist()[-1] * matrix.matrix[condition]['DW'].dropna().tolist()[-1]
                experimental_data_carbohydrate[condition] = matrix.matrix[condition]['Carbohydrate'].dropna().tolist()[-1]

    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_carotene.keys()), y=list(experimental_data_carotene.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_carotene.keys()), y=list(experimental_data_carotene.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_carotene.keys()), y=list(data_carotene.values()), zorder=2)
    plt.ylabel(r"$\beta$-Carotene (g/gDW)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/carotene_in_house.png")

    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_carotene_concentration.keys()), y=list(experimental_data_carotene_concentration.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_carotene_concentration.keys()), y=list(experimental_data_carotene_concentration.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_carotene_conc.keys()), y=list(data_carotene_conc.values()), zorder=2)
    plt.ylabel(r"$\beta$-Carotene (g/L)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/carotene_concentration_in_house.png")

    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_lutein.keys()), y=list(experimental_data_lutein.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_lutein.keys()), y=list(experimental_data_lutein.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_lutein.keys()), y=list(data_lutein.values()), zorder=2)
    plt.ylabel(r"Lutein (g/gDW)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/lutein_in_house.png")

    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_chl.keys()), y=list(experimental_data_chl.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_chl.keys()), y=list(experimental_data_chl.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_chl.keys()), y=list(data_chl.values()), zorder=2)
    plt.ylabel(r"Chlorophyll (g/L)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/chlorophyll_in_house.png")

    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_protein.keys()), y=list(experimental_data_protein.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_protein.keys()), y=list(experimental_data_protein.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_protein.keys()), y=list(data_protein.values()), zorder=2)
    plt.ylabel(r"Protein (g/gDW)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/protein_in_house.png")

    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_lipid.keys()), y=list(experimental_data_lipid.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_lipid.keys()), y=list(experimental_data_lipid.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_lipid.keys()), y=list(data_lipid.values()), zorder=2)
    plt.ylabel(r"Lipid (g/gDW)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/lipid_in_house.png")

    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_lipid.keys()), y=list(experimental_data_lipid.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_lipid.keys()), y=list(experimental_data_lipid.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_lipid_conc.keys()), y=list(data_lipid_conc.values()), zorder=2)
    plt.ylabel(r"Lipid (g/L)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/lipid_conc_in_house.png")

    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data_carbohydrate.keys()), y=list(experimental_data_carbohydrate.values()), ci='sd')
    plt.errorbar(x=list(experimental_data_carbohydrate.keys()), y=list(experimental_data_carbohydrate.values()), yerr=error, fmt='none', color='black', capsize=4)
    ax.scatter(x=list(data_carbohydrate.keys()), y=list(data_carbohydrate.values()), zorder=2)
    plt.ylabel(r"Carbohydrate (g/gDW)")
    plt.xlabel(r"Trial")
    plt.savefig(f"{DATA_PATH}/dfba/macros/carbohydrate_in_house.png")
    print(f"Total error from set of parameters: {total_error}")




def run_condition():
    """
    Runs a single condition.
    Returns
    -------

    """
    initial_parameters = json.load(open(f"{DATA_PATH}/dfba/inputs/initial_parameters.json", "r"))
    # initial_parameters = json.load(open(f"{DATA_PATH}/dfba/optimized_parameters.json", "r"))
    # model = read_model()
    # matrix = ExpMatrix(f"{DATA_PATH}/experimental/Matriz- DCCR Dunaliella salina_dfba.xlsx")
    # matrix.conditions = "Resume"
    create_dfba_model("7", initial_parameters, create_plots=True)


def run_all():
    """
    Runs all the conditions.
    Returns
    -------

    """
    initial_parameters = json.load(open(f"{DATA_PATH}/dfba/inputs/initial_parameters.json", "r"))

    # model = read_model()
    # matrix = ExpMatrix(f"{DATA_PATH}/experimental/Matriz- DCCR Dunaliella salina_dfba.xlsx")
    # matrix.conditions = "Resume"

    conditions_names = tuple(set(matrix.matrix.keys()) - {"Resume"})

    with tqdm(total=len(conditions_names), desc="Running initial conditions") as pbar:
        for condition in conditions_names:
            pbar.set_description(f"Running {condition}")
            create_dfba_model(condition, initial_parameters, create_plots=True)
            pbar.update(1)


if __name__ == '__main__':
    parameter_optimization()
    # optimize_simple_parameters()
    # run_all()
    # run_all_parallel()
    ## create random sample of conditions
