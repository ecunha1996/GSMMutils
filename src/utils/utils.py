import numpy as np


def get_productivity(data, time=48):
    previous_dw = 0
    productivity_g_l = [0]
    for i, row in data.iterrows():
        if i == 0:
            previous_dw = row["DW"]
        if i > 0:
            productivity_g_l.append((row["DW"] - previous_dw)/time)
            previous_dw = row["DW"]
    return productivity_g_l


def get_fixation(data, m_co2, m_c, carbon_in_biomass):
    fixation = []
    for i, row in data.iterrows():
        fixation.append(row["Productivity (g/L.h)"] * carbon_in_biomass * (m_co2 / m_c))
    return fixation


def get_uptake(data, m_co2):
    previous_dw = 0
    uptake = [0]
    for i, row in data.iterrows():
        if i == 0:
            previous_dw = row["DW"]
        if i > 0:
            avv = (row["DW"] + previous_dw)/2
            uptake.append((row["CO2 Fixation (gCO2/L.h)"] / avv) * 1000 / m_co2)
            previous_dw = row["DW"]
    return uptake


def get_average_uptake(data, exponential_phase) -> float:
    time_points = np.arange(exponential_phase[0]+1, exponential_phase[1]+1)
    sub_data = data.loc[data['Time (d)'].isin(time_points)]
    return np.mean(sub_data["CO2 uptake mmolCO2/gDW.h"]).round(3)