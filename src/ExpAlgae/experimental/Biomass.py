from ExpAlgae.io.reader import read_excel
from ExpAlgae.utils.utils import normalize, convert_mmol_mol_to_g_molMM


class Biomass:
    def __init__(self, biomass_reaction, filename=None):
        self.biomass_reaction = biomass_reaction
        self._standard_biomass = {}
        self.biomass_matrix = None
        if filename:
            self.load(filename)

    @property
    def standard_biomass(self):
        return self._standard_biomass

    @standard_biomass.setter
    def standard_biomass(self, biomass_composition):
        self._standard_biomass = biomass_composition

    def load(self, filename):
        self.biomass_matrix = read_excel(filename, index_col=0, sheet_name=None)

