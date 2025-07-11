import pickle
from typing import Union

import pandas as pd

from ..io import (read_matrix)
from ..io import write_matrix
from ..utils.utils import get_molecular_weight, get_element_in_biomass, get_productivity, get_fixation, \
    get_uptake, get_average_uptake, get_maximum_productivity, get_growth_rate


class ExpMatrix:
    def __init__(self, experimental_filename, model=None, conditions: str = None):
        self._conditions = None
        self._substrate_uptake_hours, self._substrate_uptake_days = {}, {}
        self.experimental_filename = experimental_filename
        self.matrix = None
        self.model = model
        self.exponential_phases = {}
        self.load()
        if conditions:
            self.conditions = conditions

    @property
    def substrate_uptake_days(self):
        return self._substrate_uptake_days

    @substrate_uptake_days.setter
    def substrate_uptake_days(self, value):
        self._substrate_uptake_days = value

    @property
    def substrate_uptake_hours(self):
        return self._substrate_uptake_hours

    @substrate_uptake_hours.setter
    def substrate_uptake_hours(self, value):
        self._substrate_uptake_hours.update(value)
        for key, value in self._substrate_uptake_hours.items():
            self._substrate_uptake_days[key] = {}
            for trial_name, uptake in value.items():
                self._substrate_uptake_days[key][trial_name] = round(uptake * 24, 4)
        new_df = pd.DataFrame.from_dict(
            data={key: value for key, value in self._substrate_uptake_days.items()})
        self.conditions = pd.concat([self.conditions.drop([column for column in new_df.columns if column in self.conditions.columns], axis=1), new_df], axis=1).rename_axis('Trial')

    def load(self):
        self.matrix = read_matrix(self.experimental_filename, sheet_name=None, index_col=0, engine="openpyxl")

    @property
    def conditions(self) -> pd.DataFrame:
        return self._conditions

    @conditions.setter
    def conditions(self, value: Union[pd.DataFrame, str]):
        if type(value) is not pd.DataFrame:
            value = self.matrix[value]
        self._conditions = value

    @conditions.getter
    def conditions(self) -> pd.DataFrame:
        return self._conditions

    def remove_trials(self, trials):
        for trial in trials:
            if trial in self.matrix.keys():
                self.matrix.pop(trial)
        self._conditions = self._conditions.drop([trial for trial in trials if trial in self._conditions.index])

    def set_exponential_phases(self, exponential_phases):
        self.exponential_phases = exponential_phases

    def save(self, filename=None):
        if filename is None:
            filename = self.experimental_filename.replace(".xlsx", "_new.xlsx")
        to_write = self.matrix
        to_write["Resume"] = self.conditions
        write_matrix(to_write, filename)
        pickle.dump(self, open(filename.replace(".xlsx", ".pkl"), "wb"))

    def get_substrate_uptake_from_biomass(self, element, substrate, header=None):
        reactions_ids = set([reaction.id for reaction in self.model.reactions])
        if not header:
            header = substrate
        m_substrate, m_element = get_molecular_weight(substrate), get_molecular_weight(element)
        temp_dict = {header: {}}
        for trial_name, data in self.matrix.items():
            if not trial_name.startswith("fachet") and not trial_name.startswith("Xi") and not trial_name.startswith("Yimei"):
                if {f"e_Biomass_trial{trial_name}__cytop"}.issubset(reactions_ids):
                    biomass_reaction_id = f"e_Biomass_trial{trial_name}__cytop"
                else:
                    biomass_reaction_id = "e_Biomass__cytop"
                carbon_in_biomass = get_element_in_biomass(self.model, element, biomass_reaction_id)
                temp_dict[header][trial_name] = self.get_substrate_uptake_for_trial(substrate, trial_name, data,
                                                                                    m_substrate, m_element,
                                                                                    carbon_in_biomass)
        self.substrate_uptake_hours = temp_dict

    def get_substrate_uptake_for_trial(self, substrate, trial_name, data, m_substrate, m_element, carbon_in_biomass):
        if 'Productivity (g/L.h)' not in data.columns:
            data['Productivity (g/L.h)'] = get_productivity(data, 48)
        data[f'{substrate} Fixation (g{substrate}/L.h)'] = get_fixation(data, m_substrate, m_element, carbon_in_biomass)
        data[f'{substrate} uptake (mmol{substrate}/gDW.h)'] = get_uptake(data, m_substrate, substrate)
        return get_average_uptake(data, self.exponential_phases[trial_name], substrate)

    def get_substrate_uptake(self, substrate, header=None):
        if not header:
            header = f"{substrate} uptake (mmol{substrate}/gDW.h)"
        temp_dict = {}
        for trial_name, data in self.matrix.items():
            if not trial_name.startswith("fachet") and not trial_name.startswith("Xi") and not trial_name.startswith("Yimei"):
                init_concentration = self.conditions.loc[trial_name, substrate]
                x_u_init = (self.matrix[trial_name].DW[str(self.exponential_phases[trial_name][0])] /
                            self.conditions['growth_rate'][trial_name])
                x_u_final = (self.matrix[trial_name].DW[str(self.exponential_phases[trial_name][1])] /
                             self.conditions['growth_rate'][trial_name])
                temp_dict[trial_name] = init_concentration / (x_u_final - x_u_init)
        self.conditions[header] = pd.Series(temp_dict)
        return temp_dict

    def get_experimental_data(self, parameter):
        if parameter == "growth_rate":
            return self.get_growth_rate()
        elif parameter == "productivity":
            return self.get_maximum_productivity()
        elif parameter == "biomass":
            return self.get_biomass()
        elif parameter == "all":
            return self.get_growth_rate(), self.get_maximum_productivity(), self.get_biomass()
        else:
            raise ValueError(f"Parameter {parameter} not recognized")

    def get_growth_rate(self):
        temp_dict = {}
        for trial_name, data in self.matrix.items():
            try:
                temp_dict[trial_name] = get_growth_rate(data, self.exponential_phases[trial_name])
            except KeyError:
                temp_dict[trial_name] = 0
        self.conditions["growth_rate"] = pd.Series(temp_dict)
        return temp_dict

    def get_maximum_productivity(self):
        temp_dict = {}
        for trial_name, data in self.matrix.items():
            try:
                temp_dict[trial_name] = get_maximum_productivity(data, self.exponential_phases[trial_name])
            except KeyError:
                temp_dict[trial_name] = 0
        self.conditions["Productivity (g/L.h)"] = pd.Series(temp_dict)
        return temp_dict

    def get_biomass(self):
        temp_dict = {}
        for trial_name, data in self.matrix.items():
            try:
                temp_dict[trial_name] = data.DW.iloc[-1]
            except KeyError:
                temp_dict[trial_name] = 0
        self.conditions["Biomass (gDW/L)"] = pd.Series(temp_dict)
        return temp_dict
