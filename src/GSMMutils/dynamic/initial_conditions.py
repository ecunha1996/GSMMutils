



def get_initial_conditions(matrix, conditions):
    active_biomass_percent = 1- (0.05 + 0.15 + 0.0131 + matrix.matrix["Resume"]["caro0"].loc[conditions] + matrix.matrix["Resume"]["chl0"].loc[conditions])
    initial_conditions = {
        "Biomass": matrix.matrix[conditions]["DW"][0],
        "ActiveBiomass": matrix.matrix[conditions]["DW"][0] * active_biomass_percent,
        "Nitrate": matrix.conditions["[N] mmol"].loc[conditions],
        "Phosphate": matrix.conditions["[P] mmol"].loc[conditions],
        "Starch": 0.05,
        "Nitrogen_quota": 5.15,
        "Phosphate_quota": 0.24,
        "Carotene": matrix.matrix["Resume"]["caro0"].loc[conditions],
        "Lutein": matrix.matrix["Resume"]["lut0"].loc[conditions],
        "Chlorophyll": matrix.matrix["Resume"]["chl0"].loc[conditions],
        "Glycerol": 0.15,
        "TAG": 0.0131
    }
    return initial_conditions