from ase.atoms import Atoms

class Descriptor:
    
    def __init__(self):
        pass
    
    def to_tensor(self, atoms: Atoms):
        raise NotImplementedError("Descriptor to tensor method is not implemented.")

    def to_derivative_tensor(self, atoms: Atoms):
        raise NotImplementedError("Descritor to derivative tensor method is not implemeneted.")