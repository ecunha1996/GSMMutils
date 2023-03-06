import os
import pickle
from os.path import join

from cobra import Metabolite

from ExpAlgae import DATA_PATH
from ExpAlgae.experimental.ExpMatrix import ExpMatrix
from ExpAlgae.model.COBRAmodel import MyModel
from ExpAlgae.utils.utils import get_biomass_mass, get_element_in_biomass, get_molecular_weight


def main():
    model = MyModel("models/model_with_trials.xml", "e_Biomass__cytop")
    model.add_medium("media.xlsx", "base_medium")
    model.setup_condition("default")
    print(model.optimize().objective_value)
    carbon_in_biomass = get_element_in_biomass(model, "C", f"e_Biomass__cytop")
    with model:
        carbohydrate_copy = model.reactions.e_Carbohydrate__cytop.copy()
        carbohydrate_copy.id = 'e_Carbohydrate_no_starch__cytop'
        biomass_copy = model.reactions.e_Biomass__cytop.copy()
        biomass_copy.id = 'e_ActiveBiomass__cytop'
        model.add_reactions([carbohydrate_copy, biomass_copy])
        model.objective = 'e_ActiveBiomass__cytop'
        model.reactions.e_Carbohydrate__cytop.bounds = (0, 0)
        # remove starch from e_carbohydrate
        carbohydrate_copy.add_metabolites({model.metabolites.C00369__chlo: -carbohydrate_copy.metabolites[model.metabolites.C00369__chlo]})
        # adjust other carbohydrates

        carbohydrate_copy.add_metabolites({model.metabolites.C00052__cytop: -carbohydrate_copy.metabolites[model.metabolites.C00052__cytop]})
        carbohydrate_copy.add_metabolites({model.metabolites.C00052__cytop: -3.35762291})
        carbohydrate_copy.add_metabolites({model.metabolites.C00096__cytop: -carbohydrate_copy.metabolites[model.metabolites.C00096__cytop]})
        carbohydrate_copy.add_metabolites({model.metabolites.C00096__cytop: -1.58181346})
        carbohydrate_copy.add_metabolites({model.metabolites.C00935__cytop: -carbohydrate_copy.metabolites[model.metabolites.C00935__cytop]})
        carbohydrate_copy.add_metabolites({model.metabolites.C00935__cytop: -2.524932989})

        carbohydrate_copy.add_metabolites({model.metabolites.C00015__cytop: -carbohydrate_copy.metabolites[model.metabolites.C00015__cytop]})
        carbohydrate_copy.add_metabolites({model.metabolites.C00015__cytop: 5.8826})
        carbohydrate_copy.add_metabolites({model.metabolites.C00035__cytop: -carbohydrate_copy.metabolites[model.metabolites.C00035__cytop]})
        carbohydrate_copy.add_metabolites({model.metabolites.C00035__cytop: 1.5818})
        carbohydrate_copy.add_metabolites({model.metabolites.C00001__cytop: -carbohydrate_copy.metabolites[model.metabolites.C00001__cytop]})
        carbohydrate_copy.add_metabolites({model.metabolites.C00001__cytop: 7.4644})

        # # adjust biomass
        biomass_copy.add_metabolites({model.metabolites.e_Carbohydrate__cytop: -biomass_copy.metabolites[model.metabolites.e_Carbohydrate__cytop]})
        biomass_copy.add_metabolites({model.metabolites.e_Carbohydrate__cytop: -0.1978})
        carbon_in_biomass = get_element_in_biomass(model, "C", f"e_ActiveBiomass__cytop")
        print(model.optimize().objective_value)
        print(model.optimize().objective_value * 0.7028)
        total = abs(sum([biomass_copy.metabolites[met] for met in biomass_copy.reactants if biomass_copy.metabolites[met] > -10]))
        print(total)
        for  met in biomass_copy.reactants:
            if biomass_copy.metabolites[met] > -10:
                new_st = biomass_copy.metabolites[met]/total
                biomass_copy.add_metabolites({model.metabolites.get_by_id(met.id):  new_st - biomass_copy.metabolites[met]})
        print(sum([biomass_copy.metabolites[met] for met in biomass_copy.reactants if biomass_copy.metabolites[met] > -10]))
        print(model.optimize().objective_value)
        carbon_in_biomass = get_element_in_biomass(model, "C", f"e_ActiveBiomass__cytop")
        print(carbon_in_biomass)
        expmatrix = pickle.load(open("experimental/Matriz- DCCR Dunaliella salina_new.pkl", "rb"))
        r = expmatrix.get_substrate_uptake_for_trial("C", "23", expmatrix.matrix["23"], get_molecular_weight("CO2"), get_molecular_weight("C"), carbon_in_biomass) * 24
        print(r)
        model.exchanges.EX_C00011__dra.bounds = (-r, 1000)
        print(model.optimize().objective_value)
        model.exchanges.EX_C00011__dra.bounds = (-8.3928, 1000)
        print(model.optimize().objective_value)
        model.write("models/model_dfba.xml")




if __name__ == "__main__":
    os.chdir(DATA_PATH)
    main()