import unittest

from GSMMutils.model.model_validator import ModelValidator
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


if __name__ == '__main__':
    unittest.main()
