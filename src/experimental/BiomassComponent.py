

class BiomassComponent:
    def __init__(self, metabolite_id, stoichiometry, parent):
        self.meta_id = metabolite_id
        self.stoichiometry = stoichiometry
        self.parent = parent