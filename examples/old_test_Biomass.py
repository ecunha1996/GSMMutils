import os
from os.path import join

import cobra

from gsmmutils import DATA_PATH
from gsmmutils.experimental.Biomass import Biomass
from gsmmutils.model.COBRAmodel import MyModel


def read_model(data_directory, filename="model.xml"):
    model = MyModel(join(join(data_directory, "models"), filename), "e_Biomass__cytop")
    model.add_medium(join(data_directory, "media.xlsx"), "base_medium")
    model.reactions.ATPm__cytop.bounds = (2.85, 2.85)
    for reaction in model.reactions:
        if reaction.lower_bound <= -1000:
            reaction.lower_bound = -10000
        if reaction.upper_bound >= 1000:
            reaction.upper_bound = 10000
    cobra_config = cobra.Configuration()
    cobra_config.bounds = -10000, 10000
    model.write(f"models/model_with_media.xml")
    model.set_prism_reaction("PRISM_white_LED__extr")
    model.reactions.R00019__chlo.bounds = (0, 0)
    return model


def simulate_all_trials(biomass, model: MyModel):
    for index, trial in biomass['macromolecules'].iterrows():
        model.adjust_biomass({"e_Protein__cytop": -trial["Protein"] / 100, "e_Carbohydrate__cytop": -trial["Carbohydrate"] / 100, "e_Lipid__cytop": -trial["Lipid"] / 100},
                             suffix=f"trial{index}")
    model.write(f"models/model_with_biomass_trials.xml")


def adjust_precursors(biomass, model):
    new_values = biomass.biomass_matrix['pigments']
    for index, row in new_values.iterrows():
        model.adjust_precursors('e_Pigment__chlo',
                                {"C05306__chlo": row['Chlorophyll a (mean)'], "C05307__chlo": row['Chlorophyll b (mean)'], "C08601__chlo": row['Lutein (mean)'], "C02094__chlo": row['B-carotene (mean)']},
                                suffix=f"trial{index}")
    model.write(f"models/model_with_trials.xml")
    # model.setup_condition("4")
    # sol = model.maximize(value=False)
    # print(model.summary(sol))
    # print(model.metabolites.e_Pigment__chlo.summary(sol))
    # model.write(f"models/model_trial4.xml")


if __name__ == '__main__':
    os.chdir(DATA_PATH)
    model = read_model(DATA_PATH, "model.xml")
    print(model.optimize())
    biomass = Biomass("e_Biomass__cytop", "experimental/Biomass_exp_composition.xlsx")
    simulate_all_trials(biomass.biomass_matrix, model)
    model = read_model(DATA_PATH, "model_with_biomass_trials.xml")
    adjust_precursors(biomass, model)
