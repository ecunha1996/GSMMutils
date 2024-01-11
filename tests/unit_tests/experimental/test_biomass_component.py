import unittest

from gsmmutils.experimental.BiomassComponent import BiomassComponent


class TestModel(unittest.TestCase):

    def test_create_biomass_component(self):
        biomass_component = BiomassComponent('metabolite', 1, None)
        self.assertEqual(biomass_component.stoichiometry, 1)
        self.assertEqual(biomass_component.parent, None)
        self.assertEqual(biomass_component.children, [])

    def test_set_parent(self):
        biomass_component = BiomassComponent('metabolite', 1, None)
        biomass_component.parent = BiomassComponent('parent', -1, None)
        self.assertEqual(biomass_component.parent.id, 'parent')
        self.assertEqual(biomass_component.children, [])

    def test_set_children(self):
        biomass_component = BiomassComponent('metabolite', 1, None)
        biomass_component.children = [BiomassComponent('child', -1, None)]
        self.assertEqual(biomass_component.parent, None)
        self.assertEqual([e.id for e in biomass_component.children], ['child'])
