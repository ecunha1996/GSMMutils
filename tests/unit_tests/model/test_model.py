import os
import unittest
from os.path import dirname, join

from cobra.io import read_sbml_model

from gsmmutils.model.COBRAmodel import MyModel


class TestModel(unittest.TestCase):
    def setUp(self):
        data_path = join(dirname(__file__), '../../data')
        os.chdir(data_path)
        self.model = MyModel("model_toy_network.xml", "e_Biomass__in")

    def test_model(self):
        self.assertEqual(len(self.model.reactions), 26)
        self.assertEqual(len(self.model.metabolites), 26)
        self.assertEqual(len(self.model.genes), 0)
        self.assertEqual(len(self.model.compartments), 2)
        self.assertEqual(len(self.model.boundary), 6)


if __name__ == '__main__':
    unittest.main()
