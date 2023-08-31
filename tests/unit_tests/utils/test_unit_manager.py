import unittest

from GSMMutils.utils.unit_manager import UnitManager, Unit


class TestUnitManager(unittest.TestCase):

    def setUp(self):
        self.unit_manager = UnitManager()

    def test_create_unit_manager(self):
        self.assertIsNotNone(self.unit_manager)

    def test_create_unit(self):
        unit = Unit("test_unit", "tests", "description")
        self.assertIsNotNone(unit)

    def test_add_unit(self):
        unit = Unit("test_unit", "tests", "description")
        self.unit_manager.add_unit(unit)
        self.assertEqual(self.unit_manager.get_unit("test_unit"), unit)

    def test_get_unit(self):
        unit = Unit("test_unit", "tests", "description")
        self.unit_manager.add_unit(unit)
        self.assertEqual(self.unit_manager.get_unit("test_unit"), unit)

    def test_get_unit_list(self):
        unit = Unit("test_unit", "tests", "description")
        self.unit_manager.add_unit(unit)
        self.assertEqual(self.unit_manager.get_units(), {unit})

    def test_set_default_units(self):
        self.unit_manager.set_default_units()
        self.assertEqual(len(self.unit_manager.get_units()), 5)


if __name__ == '__main__':
    unittest.main()
