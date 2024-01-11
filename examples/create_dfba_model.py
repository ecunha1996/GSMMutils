import os
import pickle
from cobra.flux_analysis import find_blocked_reactions
from gsmmutils import DATA_PATH
from gsmmutils.experimental.ExpMatrix import ExpMatrix
from gsmmutils.model.COBRAmodel import MyModel
from gsmmutils.utils.utils import get_element_in_biomass, get_molecular_weight


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
    pigments_copy.id = 'e_Pigments_no_car_chl__cytop'
    model.add_reactions([pigments_copy])
    model.reactions.e_Pigment__chlo.bounds = (0, 0)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05306__chlo", 0)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05307__chlo", 0)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C02094__chlo", 0)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08601__chlo", 0)
    #
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08614__chlo", -1.0633)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05433__chlo", -0.0256)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C06098__chlo", -0.0164)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08591__chlo", -0.1660)
    #
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08606__chlo", -0.0032)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05432__chlo", -0.0008)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08579__chlo", -0.0014)
    # model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C20484__chlo", -0.4527)
    #
    # model.set_stoichiometry("e_ActiveBiomass__cytop", "e_Pigment__chlo", -0.0163)

    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05306__chlo", -0.0323)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05307__chlo", -0.0109)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C02094__chlo", 0)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08601__chlo", 0)

    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08614__chlo", -1.0221)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05433__chlo", -0.0246)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C06098__chlo", -0.0158)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08591__chlo", -0.1595)

    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08606__chlo", -0.0031)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C05432__chlo", -0.0008)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C08579__chlo", -0.0014)
    model.set_stoichiometry("e_Pigments_no_car_chl__cytop", "C20484__chlo", -0.4351)

    model.set_stoichiometry("e_ActiveBiomass__cytop", "e_Pigment__chlo", -0.0163)

    return model


def correct_co2_uptake(model):
    active_biomass = model.reactions.e_ActiveBiomass__cytop
    total = abs(sum([active_biomass.metabolites[met] for met in active_biomass.reactants if active_biomass.metabolites[met] > -10]))
    print(total)
    carbon_in_biomass = get_element_in_biomass(model, "C", f"e_ActiveBiomass__cytop")
    print(carbon_in_biomass)
    # exp_matrix = pickle.load(open("experimental/Matriz- DCCR Dunaliella salina_new.pkl", "rb"))
    exp_matrix = ExpMatrix("experimental/Matriz- DCCR Dunaliella salina_new.xlsx")
    r = round(exp_matrix.get_substrate_uptake_for_trial("C", "23", exp_matrix.matrix["23"], get_molecular_weight("CO2"), get_molecular_weight("C"), carbon_in_biomass) * 24, 3)
    print(r)
    model.exchanges.EX_C00011__dra.bounds = (-r, 10000)
    print(model.optimize().objective_value)
    return model


def normalize_active_biomass(model):
    active_biomass = model.reactions.e_ActiveBiomass__cytop
    total = abs(sum([active_biomass.metabolites[met] for met in active_biomass.reactants if active_biomass.metabolites[met] > -10]))
    for met in active_biomass.reactants:
        if active_biomass.metabolites[met] > -10:
            active_biomass.metabolites[met] = model.set_stoichiometry("e_ActiveBiomass__cytop", met.id, round(active_biomass.metabolites[met] / total, 5))
    total = abs(sum([active_biomass.metabolites[met] for met in active_biomass.reactants if active_biomass.metabolites[met] > -10]))
    print(total)
    return model


def main():
    import cobra
    cobra_config = cobra.Configuration()
    cobra_config.bounds = (-10000, 10000)
    model = MyModel("models/model_with_trials.xml", "e_Biomass__cytop")
    model.add_medium("media.xlsx", "base_medium")
    model.setup_condition("default")
    for reaction in model.reactions:
        if reaction.lower_bound <= -500:
            reaction.lower_bound = -10000
        if reaction.upper_bound >= 500:
            reaction.upper_bound = 10000
    blocked = find_blocked_reactions(model)
    blocked = list(set(blocked) - {"PRISM_red_LED_674nm__extr", "PRISM_red_LED_array_653nm__extr"})
    model.remove_reactions(blocked)
    print(model.optimize().objective_value)
    with model:
        model = create_active_biomass(model)
        # model = remove_starch(model)
        # model = remove_tag(model)
        # model = remove_glycerol(model)
        model = remove_pigments(model)
        model = normalize_active_biomass(model)
        # model = correct_co2_uptake(model)
        model.write("models/model_dfba_no_caro.xml")


if __name__ == "__main__":
    os.chdir(DATA_PATH)
    main()
