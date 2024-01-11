from cobra import Metabolite

class BiomassComponent(Metabolite):
    def __init__(self, metabolite, stoichiometry, parent):
        self.stoichiometry = stoichiometry
        self._parent = parent
        self._children = []
        super().__init__(metabolite)

    @property
    def stoichiometry(self):
        return self._stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, stoichiometry):
        self._stoichiometry = stoichiometry

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent):
        self._parent = parent

    @property
    def children(self):
        return self._children

    @children.setter
    def children(self, new_children):
        self._children = new_children
        for child in self._children:
            child.parent = self
