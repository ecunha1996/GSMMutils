from cobra import Metabolite


class BiomassComponent(Metabolite):
    def __init__(self, metabolite, stoichiometry, parent):
        self.stoichiometry = stoichiometry
        self.parent = parent
        self.children = []
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
        if parent:
            parent.children.append(self)

    @property
    def children(self):
        return self._children

    @children.setter
    def children(self, children):
        self._children = children
        if children:
            children.parent = self