import os

from ExpAlgae.model.FAME2Biomass import FAME2Biomass
from test_all import read_model



if __name__ == '__main__':
    directory = r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Models"
    os.chdir(directory)
    model = read_model()
    fame2biomass = FAME2Biomass(model)
    lipid_compartments_map = {"PG": "chlo", "PC": "e_r_", "PE": "e_r_",
                              "PI": "e_r_", "MGDG": "chlo", "DGDG": "chlo",
                                "SQDG": "chlo", "DGTS": "e_r_", "DAG": "e_r_",
                                "TAG": "lip", "CL": "mito"}

    fame2biomass.run(lipid_compartments_map)
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
    fame2biomass.merge_results()
