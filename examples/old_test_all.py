from cobra import flux_analysis
from tqdm import tqdm

from gsmmutils import DATA_PATH, write_simulation, MyModel
from gsmmutils.experimental.exp_matrix import ExpMatrix
from os.path import join
from gsmmutils.graphics.plot import *
from gsmmutils.stats.stats import *


def read_model(data_directory=DATA_PATH, filename="model.xml"):
    model = MyModel(join(join(data_directory, 'models'), filename), "e_Biomass__cytop")
    model.add_medium(join(data_directory, "media.xlsx"), "base_medium")
    try:
        model.reactions.e_Biomass_ht__cytop.bounds = (0, 0)
        model.reactions.R00019__chlo.bounds = (0, 0)
        model.reactions.R00019__mito.bounds = (0, 0)
    except:
        pass
    model.set_prism_reaction("PRISM_white_LED__extr")
    model.reactions.ATPm__cytop.bounds = (1.50*24, 1.50*24)
    # blocked = flux_analysis.find_blocked_reactions(model, open_exchanges=True)
    # model.remove_reactions(blocked)
    # model.write(join(data_directory, "models/model_with_no_blocked.xml"))
    return model


def experimental_data_processing(data_directory, filename, model):
    matrix = ExpMatrix(join(data_directory, filename), model)
    matrix.conditions = "Resume"
    matrix.remove_trials(["Resume", "area", "19", "21", "R21", "25", "PC1"])
    matrix.set_exponential_phases({"1": (2, 8), "2": (2, 8), "3": (2, 16), "4": (2, 10), "5": (2, 8), "6": (2, 8), "7": (2, 16), "8": (2, 16), "9": (2, 8), "10": (2, 8), "11": (2, 12),
                                   "12": (2, 10), "13": (2, 8), "14": (2, 8), "15": (2, 14), "16": (2, 14), "17": (2, 12), "18": (2, 12), "20": (2, 14),
                                   "22": (2, 14), "23": (2, 10), "24": (2, 10), "PC1": (2, 10), "PC2": (2, 14), "PC3": (2, 10), "PC4": (2, 10), "RPC1": (2, 10), "RPC2": (2, 10),
                                   "RPC3": (2, 10), "N1": (0, 12), "N2": (0, 18), "N3": (0, 14), "N4": (0, 16), "N5": (0, 12), "N6": (0, 12), "N7": (0, 12), "N8": (0, 16), "N9": (0, 16)})
    matrix.get_experimental_data(parameter='all')
    matrix.get_substrate_uptake_from_biomass("C", "CO2", header="C00011")
    matrix.get_substrate_uptake_from_biomass("P", "C00009")
    matrix.get_substrate_uptake_from_biomass("N", "C00244")
    matrix.get_substrate_uptake("[P] mmol", header='HPO4')
    matrix.get_substrate_uptake("[N] mmol", header='NO3')
    matrix.save()
    return matrix

def simulation_for_conditions(model, conditions_df, growth_rate_df, save_in_file=False, filename=None, objective=None):
    as_dict = conditions_df.to_dict(orient='index')
    growth_rate = growth_rate_df.to_dict(orient='index')
    complete_results = {}
    error_sum = 0
    values_for_plot = {}
    model.exchanges.EX_C00011__dra.bounds = (-1000, 1000)
    for index, condition in tqdm(as_dict.items(), colour="blue"):
        if {f"e_Biomass_trial{index}__cytop"}.issubset(model.reaction_ids):
            model_copy = model.copy()
            for reaction in model_copy.reactions:
                if ("Biomass" in reaction.id and "EX_" not in reaction.id
                        and reaction.id != f"e_Biomass_trial{index}__cytop"):
                    reaction.bounds = (0, 0)
            model_copy.reactions.get_by_id(f"e_Biomass_trial{index}__cytop").bounds = (0, 1000)
            if objective:
                [setattr(x, 'objective_coefficient', 0) for x in model.reactions if x.objective_coefficient != 0]
                model_copy.reactions.get_by_id(f"e_Biomass_trial{index}__cytop").objective_coefficient = 1
                for key, value in objective.items():
                    model_copy.reactions.get_by_id(key).objective_coefficient = value
            else:
                model_copy.objective = f"e_Biomass_trial{index}__cytop"
            for met, lb in condition.items():
                lb = -lb if lb < 0 else lb
                model_copy.reactions.get_by_id("EX_" + met + "__dra").bounds = (round(-lb, 4), 1000)
            sol = model_copy.optimize()
            biomass = round(sol[f"e_Biomass_trial{index}__cytop"], 3)
            error_sum += abs(growth_rate[index]['growth_rate'] - biomass)
            complete_results[index] = sol
            values_for_plot[index] = (growth_rate[index]['growth_rate'], biomass)
    if save_in_file:
        write_simulation(complete_results, filename)
    return complete_results, values_for_plot, round(error_sum, 6)

def simulations(matrix, model):
    complete_results, values_for_plot_carbon_limited, error1 = simulation_for_conditions(model, matrix.conditions[["C00011"]], matrix.conditions[["growth_rate"]], save_in_file=True,
                                                                                         filename="../data/trial_simulations/carbon_limited")
    _, values_for_plot_p_limited, error2 = simulation_for_conditions(model, matrix.conditions[["C00009"]], matrix.conditions[["growth_rate"]], save_in_file=True,
                                                                     filename="../data/trial_simulations/p_limited")
    _, values_for_plot_carbon_and_p_limited, error3 = simulation_for_conditions(model, matrix.conditions[["C00011", "C00009"]], matrix.conditions[["growth_rate"]], save_in_file=True,
                                                                                filename="../data/trial_simulations/carbon_and_p_limited")
    _, values_for_plot_n_limited, error4 = simulation_for_conditions(model, matrix.conditions[["C00244"]], matrix.conditions[["growth_rate"]], save_in_file=True,
                                                                     filename="../data/trial_simulations/n_limited")
    _, values_for_plot_carbon_and_p_and_n_limited, error5 = simulation_for_conditions(model, matrix.conditions[["C00011", "C00009", "C00244"]], matrix.conditions[["growth_rate"]],
                                                                                      save_in_file=True, filename="../data/trial_simulations/carbon_and_p_and_n_limited")
    res = []
    print("Error1: ", error1)
    print("Error2: ", error2)
    print("Error3: ", error3)
    print("Error4: ", error4)
    print("Error5: ", error5)
    print(sum([error1, error2, error3, error4]))
    for index, element in values_for_plot_carbon_limited.items():
        new_value = values_for_plot_carbon_and_p_limited[index][1]
        new_value2 = values_for_plot_p_limited[index][1]
        new_value3 = values_for_plot_carbon_and_p_and_n_limited[index][1]
        new_value4 = values_for_plot_n_limited[index][1]
        res.append((element[0], element[1], new_value, new_value2, new_value3, new_value4))
    df = pd.DataFrame(data=res, index=complete_results.keys(), columns=["Exp", "in silico (C limited)", "in silico (C & P limited)", "in silico (P limited)", "in silico (C, P, and N limited)", "in silico (N limited)"])
    df.to_excel("../data/trial_simulations/simulation_summary.xlsx")
    barplot(df, to_show=False, path="../data/trial_simulations/simulation_summary.png")
    # model.set_photoautotrophy()
    # print(model.summary())
    # model.set_heterotrophy()
    # print(model.summary())
    # model.set_mixotrophy()
    # print(model.summary())


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


def simulations_max_carotene(matrix, model):
    complete_results, values_for_plot_carbon_and_p_and_n_limited, error = simulation_for_conditions(model, matrix.conditions[["C00011", "C00009", "C00244"]],
                                                                                                    matrix.conditions[["growth_rate"]], save_in_file=True,
                                                                                                    filename="carbon_and_p_and_n_limited", objective={"DM_C02094__chlo": 0.9})
    res = []
    print("Error1: ", error)
    for key, value in values_for_plot_carbon_and_p_and_n_limited.items():
        carotene = complete_results[key].fluxes["DM_C02094__chlo"]
        res.append((value[0], value[1], carotene))
    df = pd.DataFrame(data=res, index=complete_results.keys(), columns=["Exp", "in silico (C, P, and N limited)", "Carotene"])
    df.to_excel("../data/trial_simulations/simulation_summary_max_carotene.xlsx")
    barplot(df, to_show=False, path="../data/trial_simulations/simulation_summary_max_carotene.png")


if __name__ == '__main__':
    data_directory = r"../data"
    model = read_model(data_directory, filename="model_with_trials.xml")
    matrix = experimental_data_processing(data_directory, "experimental/Matriz- DCCR Dunaliella salina_dfba.xlsx", model)
    # matrix = ExpMatrix(join(data_directory, "Matriz- DCCR Dunaliella salina_new.xlsx"))
    matrix.conditions = "Resume"
    print("Simulating")
    simulations(matrix, model)
    simulations_max_carotene(matrix, model)
    # stats(data_directory, "Matriz- DCCR Dunaliella salina_new.xlsx")
