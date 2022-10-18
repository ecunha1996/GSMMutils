

from ExpAlgae.src.models.COBRAmodel import *
from ExpAlgae.src.ExpMatrix.ExpMatrix import *


def read_model():
    model = MyModel("model.xml", "e_Biomass__cytop")
    model.add_medium(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Env_conditions/Dsalina/media.xlsx", "base_medium")
    model.reactions.e_Biomass_ht__cytop.bounds = (0, 0)
    model.set_prism_reaction("PRISM_white_LED__extr")
    model.reactions.ATPm__cytop.bounds = (2.85, 2.85)
    print(),
    return model






if __name__ == '__main__':
    os.chdir(r"C:\Users\bisbii\OneDrive - Universidade do Minho\Models")
    model = read_model()
    matrix = ExpMatrix(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Env_conditions\Dsalina\Matriz- DCCR Dunaliella salina.xlsx", model)
    matrix.set_conditions("Resume")
    matrix.remove_trials(["Resume", "area", "19", "21", "R21"])
    matrix.set_exponential_phases({"1": (2, 8), "2": (2, 8), "3": (2, 16), "4": (2, 10), "5": (2, 8), "6": (2, 8), "7": (2,16), "8": (2, 16), "9": (2, 8),"10": (2, 8), "11": (2, 12),
                                   "12": (2, 10), "13": (2, 8) ,"14": (2, 8), "15": (2, 14), "16": (2, 14), "17": (2, 12), "18": (2, 12), "20": (2, 14),
                                   "22": (2, 14), "23": (2, 10), "24": (2, 10), "PC1": (2, 10), "PC2": (2, 14), "PC3": (2, 10), "PC4": (2, 10), "RPC1": (2, 10), "RPC2": (2, 10),
                                   "RPC3": (2, 10)})
    matrix.get_carbon_uptake()
    model.exchanges.EX_C00011__dra.lower_bound = -8.38
    print(model.optimize().objective_value)