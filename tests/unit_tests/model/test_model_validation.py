import unittest

import cobra
from pandas import DataFrame
from gsmmutils.model.model_validator import ModelValidator
from unit_tests.model.test_model import TestModel


class TestModelValidator(TestModel, unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.validator = ModelValidator(self.model, "ATPm__cytop", "__cytop")

    def test_check_balance(self):
        balance = self.validator.check_balance()
        self.assertIsInstance(balance, dict)
        self.assertEqual(balance, {})

    def test_check_energy_producing_cycles(self):
        energy_cycles = self.validator.check_energy_producing_cycles()
        self.assertIsInstance(energy_cycles, dict)
        self.assertEqual(energy_cycles, {})

    def test_check_biomass_production(self):
        biomass_production = self.validator.check_biomass_production()
        self.assertIsInstance(biomass_production, cobra.Solution)

    def test_check_blocked_reactions(self):
        blocked_reactions = self.validator.check_blocked_reactions()
        self.assertIsInstance(blocked_reactions, tuple)
        self.assertIsInstance(blocked_reactions[0], list)
        self.assertIsInstance(blocked_reactions[1], list)

    def test_check_unbounded_flux(self):
        unbounded_flux = self.validator.check_unbounded_flux()
        self.assertIsInstance(unbounded_flux, DataFrame)
        assert unbounded_flux.empty

    def test_check_sbc(self):
        sbc = self.validator.check_sbc()
        self.assertIsInstance(sbc, DataFrame)
        assert sbc.empty

    def test_check_reactions_equal_metabolites(self):
        reactions_equal_metabolites = self.validator.check_reactions_equal_metabolites()
        self.assertIsInstance(reactions_equal_metabolites, list)


if __name__ == '__main__':
    unittest.main()
