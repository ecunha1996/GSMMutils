import copy
import json
from collections import OrderedDict
from typing import List

from gsmmutils.dynamic.initial_conditions import get_initial_conditions
from gsmmutils.dynamic.rhs import get_bounds
from gsmmutils.experimental.ExpMatrix import ExpMatrix
from gsmmutils.utils.utils import get_parameter_range


class SensitivityAnalysis:
    def __init__(self, matrix: ExpMatrix = None):
        self._original_state = copy.deepcopy(self.__dict__)
        self.parameters = {}
        self.matrix = matrix

    def load_parameters(self, filename):
        with open(filename, "r") as json_file:
            self.parameters = json.load(json_file)
        self.parameters.update(get_initial_conditions(self.matrix, "1"))
        self.parameters["Eo"] = self.matrix.conditions["Light (umol/m^2.s)"].loc["1"]
        self.parameters["nacl"] = self.matrix.conditions["Salinity g/L"].loc["1"]
        self.parameters["Lr"] = self.matrix.conditions["Lr"].loc["1"]
        self.parameters["aeration"] = self.matrix.conditions["Aeration rate"].loc["1"]
        self.parameters["chlorophyll"] = 0.0063
        self.parameters = {key: Parameter(key, value, self.parameters) for key, value in self.parameters.items()}
        self.update_aliases()
        self.update_parameters()

    def update_aliases(self):
        parameters_map = OrderedDict({
            "n_quota": "Nitrogen_quota",
            "p_quota": "Phosphate_quota",
            "carotene": "Carotene",
            # "glycerol": "Glycerol",
            # "tag": "TAG",
            # "starch": "Starch",
            # "chlorophyll": "Chlorophyll",
            "X": "Biomass",
            "F": "ActiveBiomass"
        })
        for key, value in parameters_map.items():
            self.parameters[key] = Parameter(key, self.parameters[value].value, self.parameters)

    def update_parameters(self):
        ex, ex0 = get_bounds("light", {key: value.value for key, value in self.parameters.items()})
        self.parameters["Ex"], self.parameters["Ex0"] = Parameter("Ex", ex, self.parameters), Parameter("Ex0", ex0, self.parameters)
        self.parameters["q"] = Parameter("q", self.parameters["n_quota"].value / self.parameters["wNmax"].value, self.parameters)
        self.parameters["n_quota"].add_dependent(self.parameters["q"], 'self.parameters["n_quota"].value / self.parameters["wNmax"].value')
        self.parameters["n"] = Parameter("n", 1 - (self.parameters["q"].value / (self.parameters["q"].value + self.parameters["K_nitrogen_quota"].value)), self.parameters)
        self.parameters["q"].add_dependent(self.parameters["n"], '1 - (self.parameters["q"].value / (self.parameters["q"].value + self.parameters["K_nitrogen_quota"].value))')
        # self.parameters["x_storage"] = Parameter("x_storage", self.parameters["carotene"].value + self.parameters["glycerol"].value + self.parameters["tag"].value + self.parameters["starch"].value, self.parameters)
        # self.parameters["cell_size_increase"] = Parameter("cell_size_increase", 1 / (1 - self.parameters["x_storage"].value), self.parameters)
        # self.parameters["x_storage"].add_dependent(self.parameters["cell_size_increase"], '1 / (1 - self.parameters["x_storage"].value)')
        # self.parameters["z"] = Parameter("z", (self.parameters["cell_size_increase"].value - 1) / (self.parameters["t_max"].value - 1), self.parameters)
        # self.parameters["cell_size_increase"].add_dependent(self.parameters["z"], '(self.parameters["cell_size_increase"].value - 1) / (self.parameters["t_max"].value - 1)')
        self.parameters["nitrogen_mass_quota"] = Parameter("nitrogen_mass_quota", self.parameters["n_quota"].value * 14.01 / 1000, self.parameters)
        self.parameters["phosphate_mass_quota"] = Parameter("phosphate_mass_quota", self.parameters["p_quota"].value * 30.97 / 1000, self.parameters)
        self.parameters['p_quota'].add_dependent(self.parameters['phosphate_mass_quota'], 'self.parameters["p_quota"].value * 30.97 / 1000')
        self.parameters['n_quota'].add_dependent(self.parameters['nitrogen_mass_quota'], 'self.parameters["n_quota"].value * 14.01 / 1000')

    def evaluate_dynamic_expression(self, expression: str, parameter: str, param_range: List[float] = None):
        self._original_state = copy.deepcopy(self.__dict__)
        if not param_range:
            param_range = [0, self.parameters[parameter].value * 3, self.parameters[parameter].value * 3]
        param_range = get_parameter_range(param_range[0], param_range[1], param_range[2])
        res = {}
        for param in param_range:
            self.parameters[parameter].value = param
            temp_parameters = self.parameters.copy()
            for key, value in temp_parameters.items():
                if isinstance(value, Parameter):
                    temp_parameters[key] = value.value
            result = get_bounds(expression, temp_parameters)
            res[float(param)] = float(round(result, 5))
        self.reset_to_original()
        return res

    def evaluate_dynamic_expressions(self, parameters: List[str]):
        pass

    def reset_to_original(self):
        self.__dict__.update(self._original_state)


class Parameter:
    def __init__(self, name, value, parameters):
        self.name = name
        self._value = value
        self._dependent_params = {}
        self.parameters = parameters
        self.parameters[name] = self

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, new_value):
        self._value = new_value
        for dependent_param, expression in self._dependent_params.items():
            dependent_param.value = eval(expression)

    def add_dependent(self, dependent_param, expression):
        self._dependent_params[dependent_param] = expression

