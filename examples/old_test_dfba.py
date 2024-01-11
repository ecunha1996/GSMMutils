from os.path import join

import numpy as np

from gsmmutils.dynamic.soa import soa
from gsmmutils.model.COBRAmodel import MyModel
from gsmmutils.utils.utils import get_element_in_biomass, get_biomass_mass


def read_model(data_directory):
    model = MyModel(join(data_directory, "model_trial4.xml"), "e_Biomass_trial4__cytop")
    model.add_medium(join(data_directory, "media.xlsx"), "base_medium")
    model.reactions.e_Biomass_ht__cytop.bounds = (0, 0)
    model.reactions.R00019__chlo.bounds = (0, 0)
    model.reactions.R00019__mito.bounds = (0, 0)
    model.set_prism_reaction("PRISM_white_LED__extr")
    model.reactions.ATPm__cytop.bounds = (2.85, 2.85)
    return model


if __name__ == '__main__':
    data_directory = r"../data"
    model = read_model(data_directory)
    model.exchanges.EX_C00011__dra.lower_bound = -1000
    model.exchanges.EX_C00205__dra.lower_bound = -1000
    dfba = soa(model)
    wn = get_element_in_biomass(model, 'N', 'e_Biomass_trial4__cytop')
    kinetics = {"PRISM_white_LED__extr": {'type': 'light'}, 'EX_C00009__dra': {'type': 'micmen', 'Km': 0.0038, 'Vmax': 0.39*24}, 'DM_C02094__chlo': {"type": 'caro', 'n': 2, 'wn':wn}}
    metabolites_to_follow = {"C00009__extr": "EX_C00009__dra", "e_Biomass_trial4__cytop": "e_Biomass_trial4__cytop", "C00244__extr": "EX_C00244__dra", "C02094__chlo": "DM_C02094__chlo"}
    initial_conditions = {"C00009__extr": 0.15, "e_Biomass_trial4__cytop": 0.137, "C00244__extr": 12.5, "C02094__chlo": 0.0}
    dfba.set_params(np.arange(0, 16, 1), kinetics = kinetics, metabolites_to_follow=metabolites_to_follow, initial_conditions = initial_conditions)
    dfba.run()
    dfba.concentrations.T["C02094__chlo"] = dfba.concentrations.T["C02094__chlo"] / 1000 * dfba.model.metabolites.get_by_id("C02094__chlo").formula_weight / dfba.concentrations.T["e_Biomass_trial4__cytop"]
    print(dfba.concentrations)
    dfba.plot_concentrations(['e_Biomass_trial4__cytop', 'C00009__extr', 'C02094__chlo'])
