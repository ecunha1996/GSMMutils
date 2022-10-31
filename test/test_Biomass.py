from experimental.Biomass import Biomass
from test_all import read_model







if __name__ == '__main__':
    data_directory = r"../data"
    model = read_model(data_directory)
    biomass = Biomass(model, "e_Biomass__cytop")
    print(biomass.standard_biomass['e_Protein__cytop'])