from ase.atom import Atom

class AtomGraphNode():

    def __init__(self, atom: Atom):
        self.visited = False
        self.atom = atom

    def isVisited(self):
        return self.visited

    def setVisited(self, state: bool):
        self.visited = state

    def getAtom(self): return self.atom


