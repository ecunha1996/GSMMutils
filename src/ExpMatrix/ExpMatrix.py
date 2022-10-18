from cobra import Metabolite

from ExpAlgae.src.io.writer import *
from ExpAlgae.src.io.reader import read_matrix
from models.COBRAmodel import get_biomass_mass
from utils.utils import *


class ExpMatrix:
    def __init__(self, experimental_filename, model):
        self.conditions = None
        self._substrate_uptake_hours, self._substrate_uptake_days = {}, {}
        self.experimental_filename = experimental_filename
        self.matrix = None
        self.model = model
        self.exponential_phases = {}
        self.load()

    @property
    def substrate_uptake_days(self):
        return {key: round(value * 24, 3) for key, value in self.substrate_uptake_hours.items()}

    @substrate_uptake_days.setter
    def substrate_uptake_days(self, value):
        self._substrate_uptake_days = value

    @property
    def substrate_uptake_hours(self):
        return self._substrate_uptake_hours

    @substrate_uptake_hours.setter
    def substrate_uptake_hours(self, value):
        self._substrate_uptake_hours = value

    def load(self):
        self.matrix = read_matrix(self.experimental_filename, sheet_name=None)

    def set_conditions(self, conditions_sheet_name):
        self.conditions = self.matrix[conditions_sheet_name]

    def remove_trials(self, trials):
        for trial in trials:
            self.matrix.pop(trial)

    def set_exponential_phases(self, exponential_phases):
        self.exponential_phases = exponential_phases

    def save(self, filename=None):
        if filename is None:
            filename = self.experimental_filename.replace(".xlsx", "_new.xlsx")
        write_matrix(self.matrix, filename)

    def get_carbon_uptake(self):
        res = get_biomass_mass(self.model)
        percentage_map = {}
        for key in res[1]:
            percentage_map[key] = round(res[1][key] * Metabolite(formula=key).formula_weight / 1000, 3)
        carbon_in_biomass = percentage_map['C']
        m_co2 = Metabolite(formula="CO2").formula_weight
        m_c = Metabolite(formula="C").formula_weight
        for trial_name, data in self.matrix.items():
            self.get_carbon_uptake_for_trial(trial_name, data, m_co2, m_c, carbon_in_biomass)
        print(self.substrate_uptake_days)

    def get_carbon_uptake_for_trial(self, trial_name, data, m_co2, m_c, carbon_in_biomass):
        data['Productivity (g/L.h)'] = get_productivity(data, 48)
        data['CO2 Fixation (gCO2/L.h)'] = get_fixation(data, m_co2, m_c, carbon_in_biomass)
        data['CO2 uptake mmolCO2/gDW.h'] = get_uptake(data, m_co2)
        self.substrate_uptake_hours[trial_name] = get_average_uptake(data, self.exponential_phases[trial_name])
