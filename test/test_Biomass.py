from experimental.Biomass import Biomass
from test_all import read_model


def simulate_all_trials(biomass, model):
    for index, trial in biomass.biomass_matrix.iterrows():
        model.adjust_biomass({"e_Protein__cytop": -trial["Protein"]/100, "e_Carbohydrate__cytop": -trial["Carbohydrate"]/100, "e_Lipid__cytop": -trial["Lipid"]/100},
                             suffix=f"trial{index}")
        model.objective = model.reactions.get_by_id(f"e_Biomass_trial{index}__cytop")
        for reaction in model.reactions:
            if "Biomass" in reaction.id and "EX_" not in reaction.id and reaction.id != f"e_Biomass_trial{index}__cytop":
                reaction.bounds = (0, 0)
        print(model.summary())




if __name__ == '__main__':
    data_directory = r"../data"
    model = read_model(data_directory)
    biomass = Biomass("e_Biomass__cytop", "Biomass_exp_composition.xlsx")
    simulate_all_trials(biomass, model)