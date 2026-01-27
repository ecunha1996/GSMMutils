class UnitManager:
    def __init__(self):
        self._units = set()

    def add_unit(self, unit):
        self._units.add(unit)

    def remove_unit(self, unit):
        self._units.remove(unit)

    def get_unit(self, name):
        for unit in self._units:
            if unit.name == name:
                return unit
        return None

    def get_units(self):
        return self._units

    def get_units_by_type(self, unit_type):
        units = []
        for unit in self._units:
            if unit.unit_type == unit_type:
                units.append(unit)
        return units

    def set_default_units(self):
        default_unit_types = {"gDW": "mass", "gMM": "mass", "gM": "mass", "molM": "amount",
                              "molMM": "amount"}
        default_unit_descriptions = {"gDW": "gram dry weight", "gMM": "gram mass of a macromolecule",
                                     "gM": "gram mass of a precursor of a macromolecule",
                                     "mol": "mole of a precursor of a macromolecule",
                                     "molMM": "mole of a macromolecule",
                                     "molM": "mole of a precursor of a macromolecule"}

        for unit, unit_type in default_unit_types.items():
            self.add_unit(Unit(unit, unit_type, default_unit_descriptions[unit]))

    def convert(self, value, from_unit, to_unit):
        # TODO: implement conversion
        pass


class Unit:
    def __init__(self, name, unit_type, description=None):
        self._name = name
        self._unit_type = unit_type
        self._description = description

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def unit_type(self):
        return self._unit_type

    @unit_type.setter
    def unit_type(self, value):
        self._unit_type = value

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, value):
        self._description = value
