from GSMMutils.model.COBRAmodel import MyModel
from GSMMutils.model.model_validator import check_consistency
from os.path import join





if __name__ == '__main__':
    data_directory = "../data"
    model = MyModel(join(data_directory, "models/model_with_trials.xml"), "e_Biomass__cytop")
    check_consistency(model, atpm_reaction="ATPm__cytop", cytoplasm_abb="__cytop")