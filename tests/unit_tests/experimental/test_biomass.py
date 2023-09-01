import unittest
from os.path import join, dirname
from os import chdir

from GSMMutils import Biomass, MyModel


class TestBiomass(unittest.TestCase):

    def setUp(self) -> None:
        data_path = join(dirname(__file__), '../../data')
        chdir(data_path)
        self.model = MyModel("model_toy_network.xml", "e_Biomass__in")

    def test_create_biomass(self):
        biomass = Biomass(self.model.reactions.get_by_id("e_Biomass__in"), "Matrix_test.xlsx")
        self.assertEqual(biomass.biomass_reaction.id, "e_Biomass__in")
