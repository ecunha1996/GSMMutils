import os
from collections import OrderedDict
from os.path import join

from GSMMutils import DATA_PATH, MyModel
from GSMMutils.model.FAME2Biomass import FAME2Biomass, load_results

if __name__ == '__main__':
    directory = join(DATA_PATH, 'fame2biomass')
    model = MyModel(join(join(DATA_PATH, 'models'), 'model_ng.xml'), "e_Biomass__cytop")
    model.add_medium(join(DATA_PATH, "media.xlsx"), "base_medium")
    model.id = model.id.strip("_v2")
    fame2biomass = FAME2Biomass(model, data_directory=directory)
    lipid_compartments_map = {"PG": "chlo", "PC": "er", "PE": "er",
                              "PI": "er", "MGDG": "chlo", "DGDG": "chlo",
                                "SQDG": "chlo", "DGTS": "er", "DAG": "er",
                                "TAG": "lip", "CL": "mito"}
    abb_compartment_id_map = {'chlo': 'C_00005', 'er': 'C_00003', 'lip': 'C_00010', 'mito': 'C_00002'}
    # fame2biomass.batch_run(lipid_compartments_map, abb_compartment_id_map)
    # fame2biomass.run(lipid_compartments_map)
    # main(model, 'PG', compartment='chlo')
    # main(model, 'PC', compartment='e_r_')
    # main(model, 'PE', compartment='e_r_')
    # main(model, 'PI', compartment='e_r_')
    # main(model, 'MGDG', compartment='chlo')
    # main(model, 'DGDG', compartment='chlo')
    # main(model, 'SQDG', compartment='chlo')
    # main(model, 'DGTS', compartment='e_r_')
    # main(model, 'DAG', compartment='e_r_')
    # fame2biomass.run("TAG", compartment='lip', compartment_id="C_lip")
    # fame2biomass.run("PG", compartment='chlo', compartment_id="C_chlo")
    # obj_list = OrderedDict({
    #     "12_0": 0.005597227,
    #     "14_0": 0.105310116,
    #     "16_0": 0.10763085,
    #     "16_1": 0.116222095,
    #     "16_2": 0.028731998,
    #     "16_3": 0.056024689,
    #     "18_0": 0.004043437,
    #     "18_1": 0.021664232,
    #     "18_2": 0.013548449,
    #     "18_3": 0.034911874,
    #     "20_1": 0.001175418,
    #     "20_4_v2": 0.032446063,
    #     "20_5": 0.472693552,
    # })
    # fame2biomass.run("MGDG", obj_list, compartment='chlo', compartment_id="C_chlo")
    # obj_list = OrderedDict({
    #     "14_0": 0.061406026,
    #     "16_0": 0.243595112,
    #     "16_1": 0.317911444,
    #     "16_2": 0.008976943,
    #     "16_3": 0.00515868,
    #     "18_0": 0.004380012,
    #     "18_1": 0.02062296,
    #     "18_2": 0.021515751,
    #     "18_3": 0.010290946,
    #     "20_4_v2": 0.007786687,
    #     "20_5": 0.29835544,
    # })
    # fame2biomass.run("DGDG", obj_list, compartment='chlo', compartment_id="C_chlo")
    # #main(model, 'CL', compartment='mito')
    # obj_list = OrderedDict({
    #     "14_0": 0.043577258,
    #     "16_0": 0.371276242,
    #     "16_1": 0.285966522,
    #     "16_2": 0.00680635,
    #     "16_3": 0.009150848,
    #     "18_0": 0.00532333,
    #     "18_1": 0.11922693,
    #     "18_2": 0.042867273,
    #     "18_3": 0.045057112,
    #     "20_1": 0.00325559,
    #     "20_2": 0.001583801,
    #     "20_3_v2": 0.002199723,
    #     "20_4_v2": 0.0270126,
    #     "20_5": 0.03669642,
    # })
    # fame2biomass.run("SQDG", obj_list, compartment='chlo', compartment_id="C_chlo")
    # obj_list = OrderedDict({
    #     "14_0": 0.040135947,
    #     "16_0": 0.130508986,
    #     "16_1": 0.182152532,
    #     "17_0": 0.004561512,
    #     "17_1": 0.0224924,
    #     "18_0": 0.00436691,
    #     "18_1": 0.019950212,
    #     "18_2": 0.022195869,
    #     "18_3": 0.00988604,
    #     "18_4": 0.008658579,
    #     "20_4": 0.058043995,
    #     "20_5": 0.497047018,
    # })
    # fame2biomass.run("DGTS", obj_list, compartment='er', compartment_id="C_er")
    # main(model, 'CL', compartment='mito')
    # fame2biomass.merge_results()
    # load_results()

    obj_list = OrderedDict({
        "14_0": 0.025498544,
        "16_0": 0.29253801,
        "16_1": 0.391419873,
        "16_2": 0.010925711,
        "16_3": 0.002817519,
        "18_0": 0.006762045,
        "18_1": 0.064251166,
        "18_2": 0.065408175,
        "18_3": 0.033763266/2,
        "18_3v2": 0.033763266 / 2,
        "20_3": 0.031403594,
        "20_4": 0.021741853,
        "20_5": 0.053470244,
    })
    fame2biomass.run("PC", obj_list, compartment='er', compartment_id="C_er")

