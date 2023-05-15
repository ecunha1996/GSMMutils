import os
from os.path import join

from ExpGSMM import DATA_PATH, MyModel
from ExpGSMM.model.FAME2Biomass import FAME2Biomass
from test_all import read_model


if __name__ == '__main__':
    directory = join(DATA_PATH, 'fame2biomass')
    model = MyModel(join(join(DATA_PATH, 'models'), 'model.xml'), "e_Biomass__cytop")
    model.add_medium(join(DATA_PATH, "media.xlsx"), "base_medium")
    fame2biomass = FAME2Biomass(model, data_directory=directory)
    lipid_compartments_map = {"PG": "chlo", "PC": "er", "PE": "er",
                              "PI": "er", "MGDG": "chlo", "DGDG": "chlo",
                                "SQDG": "chlo", "DGTS": "er", "DAG": "er",
                                "TAG": "lip", "CL": "mito"}
    abb_compartment_id_map = {'chlo': 'C_00005', 'er': 'C_00003', 'lip': 'C_00010', 'mito': 'C_00002'}
    fame2biomass.batch_run(lipid_compartments_map, abb_compartment_id_map)
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
    # main(model, 'TAG', compartment='lip')
    #main(model, 'CL', compartment='mito')
    # fame2biomass.merge_results()
