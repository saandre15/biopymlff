from ase.atom import Atom

class AtomGraphNode():

    def __init__(self, atom: Atom):
        self.visited = False
        self.atom = atom
        self.organic_unpaired_electrons = {
            'H': 0,
            'C': 0,
            'N': 2,
            'O': 4,
            'F': 6,
            'Cl': 6,
            'Br': 6,
            'I': 6,
            'Si': 0,
            'P': 2,
            'S': 4
        }

    def isVisited(self):
        return self.visited

    def setVisited(self, state: bool):
        self.visited = state

    def getAtom(self): return self.atom


