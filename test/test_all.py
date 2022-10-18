
from src.models.COBRAmodel import *
from src.ExpMatrix.ExpMatrix import *
from os.path import join
from graphics.plot import *


def read_model(data_directory):
    model = MyModel(join(data_directory, "model.xml"), "e_Biomass__cytop")
    model.add_medium(join(data_directory, "media.xlsx"), "base_medium")
    model.reactions.e_Biomass_ht__cytop.bounds = (0, 0)
    model.set_prism_reaction("PRISM_white_LED__extr")
    model.reactions.ATPm__cytop.bounds = (2.85, 2.85)
    print(),
    return model



if __name__ == '__main__':
    data_directory = r"../data"
    model = read_model(data_directory)
    matrix = ExpMatrix(join(data_directory, "Matriz- DCCR Dunaliella salina.xlsx"), model)
    matrix.conditions = "Resume"
    matrix.remove_trials(["Resume", "area", "19", "21", "R21", "25", "PC1"])
    matrix.set_exponential_phases({"1": (2, 8), "2": (2, 8), "3": (2, 16), "4": (2, 10), "5": (2, 8), "6": (2, 8), "7": (2,16), "8": (2, 16), "9": (2, 8),"10": (2, 8), "11": (2, 12),
                                   "12": (2, 10), "13": (2, 8) ,"14": (2, 8), "15": (2, 14), "16": (2, 14), "17": (2, 12), "18": (2, 12), "20": (2, 14),
                                   "22": (2, 14), "23": (2, 10), "24": (2, 10), "PC1": (2, 10), "PC2": (2, 14), "PC3": (2, 10), "PC4": (2, 10), "RPC1": (2, 10), "RPC2": (2, 10),
                                   "RPC3": (2, 10)})
    matrix.get_experimental_data(parameter='all')
    matrix.get_substrate_uptake_from_biomass("C", "CO2", header = "C00011")
    matrix.get_substrate_uptake_from_biomass("P", "HPO4")
    matrix.get_substrate_uptake("[P] mmol", header = 'C00009')
    matrix.get_substrate_uptake("[N] mmol", header='C00244')
    matrix.save()
    complete_results, values_for_plot_carbon_limited, _ = simulation_for_conditions(model, matrix.conditions[["C00011"]], matrix.conditions[["growth_rate"]], save_in_file=True, filename="carbon_limited")
    _, values_for_plot_carbon_and_p_limited, _ = simulation_for_conditions(model, matrix.conditions[["C00011", "C00009"]], matrix.conditions[["growth_rate"]], save_in_file=True, filename="carbon_and_p_limited")
    _, values_for_plot_carbon_and_p_and_n_limited, _ = simulation_for_conditions(model, matrix.conditions[["C00011", "C00009", "C00244"]], matrix.conditions[["growth_rate"]], save_in_file=True, filename="carbon_and_p_and_n_limited")
    res = []
    for index, element in enumerate(values_for_plot_carbon_limited):
        new_value = values_for_plot_carbon_and_p_limited[index][1]
        new_value2 = values_for_plot_carbon_and_p_and_n_limited[index][1]
        res.append((element[0], element[1], new_value))
    df = pd.DataFrame(data=res, index= complete_results.keys(), columns=["Exp", "in silico (C limited)", "in silico (C & P limited)"])
    df.to_excel("simulation_summary.xlsx")
    barplot(df)