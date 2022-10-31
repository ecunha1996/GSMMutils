from os import path
from model.COBRAmodel import MyModel
from utils.utils import infer_biomass_from_model


class Biomass:
    def __init__(self, model, biomass_reaction, set_automatic_biomass=True):
        self.model = model
        self.biomass_reaction = biomass_reaction
        if set_automatic_biomass:
            self.set_Biomass()

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        if type(model) == MyModel:
            self._model = model
        elif type(model) == path or type(model) == str:
            self._model = MyModel(model)
        else:
            raise TypeError("model must be a MyModel object or a path to a model file")

    @property
    def standard_biomass(self):
        return self._standard_biomass

    @standard_biomass.setter
    def standard_biomass(self, biomass_reaction):
        if not biomass_reaction:
            self._standard_biomass = infer_biomass_from_model(self.model, self.biomass_reaction)

    def set_Biomass(self, biomass_reaction=None):
        self.standard_biomass = biomass_reaction


