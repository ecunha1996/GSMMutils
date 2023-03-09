import os
import pickle
from os.path import join

from cobra import Metabolite
from cobra.flux_analysis import find_blocked_reactions

from ExpAlgae import DATA_PATH
from ExpAlgae.experimental.ExpMatrix import ExpMatrix
from ExpAlgae.model.COBRAmodel import MyModel
from ExpAlgae.utils.utils import get_biomass_mass, get_element_in_biomass, get_molecular_weight


def create_active_biomass(model):
    biomass_copy = model.reactions.e_Biomass__cytop.copy()
    biomass_copy.id = 'e_ActiveBiomass__cytop'
    model.add_reactions([biomass_copy])
    model.objective = 'e_ActiveBiomass__cytop'
    return model


def remove_starch(model):
    carbohydrate_copy = model.reactions.e_Carbohydrate__cytop.copy()
    carbohydrate_copy.id = 'e_Carbohydrate_no_starch__cytop'
    model.add_reactions([carbohydrate_copy])

    model.reactions.e_Carbohydrate__cytop.bounds = (0, 0)
    # remove starch from e_carbohydrate
    carbohydrate_copy.add_metabolites({model.metabolites.C00369__chlo: -carbohydrate_copy.metabolites[model.metabolites.C00369__chlo]})
    # adjust other carbohydrates
    model.set_stoichiometry("e_Carbohydrate_no_starch__cytop", "C00052__cytop", -3.35762291)
    model.set_stoichiometry("e_Carbohydrate_no_starch__cytop", "C00096__cytop", -1.58181346)
    model.set_stoichiometry("e_Carbohydrate_no_starch__cytop", "C00935__cytop", -2.524932989)
    model.set_stoichiometry("e_Carbohydrate_no_starch__cytop", "C00015__cytop", 5.8826)
    model.set_stoichiometry("e_Carbohydrate_no_starch__cytop", "C00035__cytop", 1.5818)
    model.set_stoichiometry("e_Carbohydrate_no_starch__cytop", "C00001__cytop", 7.4644)

    # # adjust biomass
    model.set_stoichiometry("e_ActiveBiomass__cytop", "e_Carbohydrate__cytop", -0.1978)
    print(model.optimize().objective_value)
    print(model.optimize().objective_value * 0.7028)
    return model


def remove_tag(model):
    lipid_copy = model.reactions.e_Lipid__cytop.copy()
    lipid_copy.id = 'e_Lipid_no_tag__cytop'
    model.add_reactions([lipid_copy])
    model.reactions.e_Lipid__cytop.bounds = (0, 0)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C00422__lip", 0)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C00641__er", -0.1607)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C00162__cytop", -0.1307)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C13508__chlo", -0.1292)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C00344__chlo", -0.2049)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C00157__er", -0.0438)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C00350__er", -0.0489)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C01194__er", -0.0373)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C03692__chlo", -0.2690)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C06037__chlo", -0.2165)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C18169__er", -0.0978)
    model.set_stoichiometry("e_Lipid_no_tag__cytop", "C05980__mito", -0.0299)
    # # adjust biomass
    model.set_stoichiometry("e_ActiveBiomass__cytop", "e_Lipid__cytop", -0.0979)
    print(model.optimize().objective_value)
    print(model.optimize().objective_value * 0.6898)
    return model


def remove_glycerol(model):
    acids_copy = model.reactions.e_Acid__cytop.copy()
    acids_copy.id = 'e_Acid_no_glycerol__cytop'
    model.add_reactions([acids_copy])
    model.reactions.e_Acid__cytop.bounds = (0, 0)
    model.set_stoichiometry("e_Acid_no_glycerol__cytop", "C00116__cytop", 0)
    model.set_stoichiometry("e_Acid_no_glycerol__cytop", "C00033__cytop", -5.5505)
    model.set_stoichiometry("e_Acid_no_glycerol__cytop", "C00246__cytop", -4.4998)
    model.set_stoichiometry("e_Acid_no_glycerol__cytop", "C00163__cytop", -3.7835)
    model.set_stoichiometry("e_ActiveBiomass__cytop", "e_Acid__cytop", -0.0067)
    return model


def remove_pigments(model):
    pigments_copy = model.reactions.e_Pigment__chlo.copy()
    pigments_copy.id = 'e_Pigments_no_car_chl_lut__cytop'
    model.add_reactions([pigments_copy])
    model.reactions.e_Pigment__chlo.bounds = (0, 0)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C05306__chlo", 0)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C05307__chlo", 0)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C02094__chlo", 0)

    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C08614__chlo", -1.0095)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C05433__chlo", -0.0243)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C06098__chlo", -0.0156)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C08591__chlo", -0.1576)

    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C08606__chlo", -0.0030)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C05432__chlo", -0.0008)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C08579__chlo", -0.0014)
    model.set_stoichiometry("e_Pigments_no_car_chl_lut__cytop", "C20484__chlo", -0.4298)

    model.set_stoichiometry("e_ActiveBiomass__cytop", "e_Acid__cytop", -0.0172)

    return model


def correct_co2_uptake(model, active_biomass):
    total = abs(sum([active_biomass.metabolites[met] for met in active_biomass.reactants if active_biomass.metabolites[met] > -10]))
    print(total)
    # for met in active_biomass.reactants:
    #     if active_biomass.metabolites[met] > -10:
    #         new_st = active_biomass.metabolites[met] / total
    #         active_biomass.add_metabolites({model.metabolites.get_by_id(met.id): new_st - active_biomass.metabolites[met]})
    # print(round(abs(sum([active_biomass.metabolites[met] for met in active_biomass.reactants if active_biomass.metabolites[met] > -10])), 4))
    # print(model.optimize().objective_value)
    carbon_in_biomass = get_element_in_biomass(model, "C", f"e_ActiveBiomass__cytop")
    print(carbon_in_biomass)
    expmatrix = pickle.load(open("experimental/Matriz- DCCR Dunaliella salina_new.pkl", "rb"))
    r = round(expmatrix.get_substrate_uptake_for_trial("C", "23", expmatrix.matrix["23"], get_molecular_weight("CO2"), get_molecular_weight("C"), carbon_in_biomass) * 24, 3)
    print(r)
    model.exchanges.EX_C00011__dra.bounds = (-r, 1000)
    print(model.optimize().objective_value)
    return model


def main():
    model = MyModel("models/model_with_trials.xml", "e_Biomass__cytop")
    model.add_medium("media.xlsx", "base_medium")
    model.setup_condition("default")
    blocked = find_blocked_reactions(model)
    model.remove_reactions(blocked)
    print(model.optimize().objective_value)
    with model:
        model = create_active_biomass(model)
        model = remove_starch(model)
        model = remove_tag(model)
        model = remove_glycerol(model)
        model = remove_pigments(model)
        model = correct_co2_uptake(model, model.reactions.e_ActiveBiomass__cytop)
        model.write("models/model_dfba.xml")




if __name__ == "__main__":
    os.chdir(DATA_PATH)
    main()