import pandas as pd

from ..graphics.plot import plot_two_axis
from ..utils.utils import get_light_kinetics, get_micmen_kinetics, get_caro_kinetics


class soa:
    def __init__(self, model):
        self.fluxes = None
        self.concentrations = None
        self.metabolites_to_follow_reverse = None
        self.metabolites_to_follow = None
        self.time_span = None
        self._time_span = None
        self._kinetics = {}
        self.model = model
        self.initial_conditions = None
        self.parameters = {}

    @property
    def timestamps(self):
        return self._time_span

    @timestamps.setter
    def timestamps(self, time_span):
        self._time_span = time_span

    @property
    def kinetics(self):
        return self._kinetics

    @kinetics.setter
    def kinetics(self, kinetics):
        self.parameters = kinetics
        for key, value in kinetics.items():
            if key in self.model.reactions:
                if value['type'] == 'light':
                    self.assing_light_kinetics(key)
                elif value['type'] == 'micmen':
                    self.assign_micmen_kinetics(key)
                elif value['type'] == 'caro':
                    self.assin_caro_kinetics(key)
            else:
                raise Exception("Reaction not found in model")

    def set_params(self, time_span, kinetics, metabolites_to_follow, initial_conditions):
        self.time_span = time_span
        self.kinetics = kinetics
        self.metabolites_to_follow = metabolites_to_follow
        self.metabolites_to_follow_reverse = {v: k for k, v in metabolites_to_follow.items()}
        self.initial_conditions = initial_conditions

    def update_kinetics(self, biomass, concentrations):
        for key, value in self.kinetics.items():
            if key in self.metabolites_to_follow_reverse and self.metabolites_to_follow_reverse[key] in concentrations.index:
                conc = concentrations.loc[self.metabolites_to_follow_reverse[key]].values[-1]
            else:
                conc = 0
            val = value(biomass, conc, self.parameters[key])
            self.model.reactions.get_by_id(key).bounds = val

    def run(self):
        concentrations, fluxes = pd.DataFrame.from_dict(self.initial_conditions, orient='index', columns=['0']), pd.DataFrame(index=[r.id for r in self.model.reactions])
        self.update_kinetics(self.initial_conditions[self.model.bio_reaction.id], concentrations)
        X = self.initial_conditions[self.model.bio_reaction.id]
        [setattr(x, 'objective_coefficient', 0) for x in self.model.reactions if x.objective_coefficient != 0]
        reaction_1 = self.model.reactions.get_by_id(self.model.bio_reaction.id)
        reaction_2 = self.model.reactions.DM_C02094__chlo
        reaction_1.objective_coefficient = 1
        reaction_2.objective_coefficient = -1
        for i in self.time_span:
            # solution = pfba(self.model)
            solution = self.model.optimize()
            fluxes[str(i)] = solution.fluxes
            X += solution.fluxes.loc[self.model.bio_reaction.id] * X
            car_mg = concentrations.T["C02094__chlo"].loc[str(i)] / 1000 * self.model.metabolites.get_by_id("C02094__chlo").formula_weight / X / 1000
            X += car_mg
            for key, value in self.metabolites_to_follow.items():
                if "biomass" not in key.lower():
                    new_conc = concentrations.loc[key, str(i)] + solution.fluxes.loc[value] * X
                    concentrations.loc[key, str(i + 1)] = max(new_conc, 0)
            concentrations.loc[self.model.bio_reaction.id, str(i + 1)] = X
            self.update_kinetics(X, concentrations)
        self.concentrations = concentrations
        self.fluxes = fluxes

    def assing_light_kinetics(self, key):
        self.kinetics[key] = get_light_kinetics

    def assign_micmen_kinetics(self, key):
        self.kinetics[key] = get_micmen_kinetics

    def assin_caro_kinetics(self, key):
        self.kinetics[key] = get_caro_kinetics

    def plot_concentrations(self, columns=None):
        if columns:
            plot_two_axis(self.concentrations.T[columns], secondary=['e_Biomass_trial4__cytop'])
