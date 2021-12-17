import os as os

from enum import Enum

from ase.atoms import Atoms
from ase.io import write

from pysmiles import write_smiles

from ..util.getenv import getenv
from ..data.atom_graph import AtomGraph

general_params = getenv()["general"]

class Library:
    """ 
    Data Library used by the other classes to access atomic data.
    """

    def __init__(self, path: str=None):
        self.lib_path = ""
        self.lib_cache = self.load_lib()
        self.file_count = len(self.lib_cache)

    def get_path(self):
        if not self.lib_path: raise NotImplementedError("Cannot get path for an abstract Library class.")
        return self.lib_path

    def add_to_lib(self, atoms: Atoms):
        """ 

        """        
        mol_graph = AtomGraph(atoms)
        smiles = write_smiles(molecule=mol_graph.get_graph())
        self.lib_cache[smiles] = atoms
        if not os.path.exists(self.lib_path): os.mkdir(self.lib_path)
        filename = str(atoms.symbols()) + "_" + self.file_count + ".pdb"
        write(filename, atoms) 
        self.file_count+=1

    def load_lib(self):
        files = [file for file in os.listdir(self.lib_path) if os.path.isfile(self.lib_path + "/" + file)]
        atoms = list(np.zeros(shape=(len(files), 1)))
        atoms = dict()
        self.file_count = 0
        for file in files:
            atom = read_proteindatabank(file)
            atoms[counter] = atom
            self.file_count+=1
        return atoms
    
    def in_lib(self, atoms: Atoms):
        mol_graph = AtomGraph(atoms)
        smiles = write_smiles(molecule=mol_graph)
        mol: Atoms = self.lib_cache[smiles]
        return mol == atoms

    def extract_atoms(dataset_size: int) -> list: # NOTE: These are by subsystems size not fragment size
        """ 
        Randomly extract previously added atoms based off subsystem size
        """
        subsystems = []
        # TODO: Loop through all of the atoms and get lsqc to fragment and create the subsystems. 
        # Purposely allow it to fail and read the subsystems and then add it to any array until the dataset size requirement is met.

        raise Exception("Library size is not large enough. Add more atoms to the library before calling this method.")

    

class ProteinLibrary(Library):
    def get_path(self):
        return general_params["protein_lib_path"] if not self.lib_path else self.lib_path


class DnaRnaLibrary(Library):
    def get_path(self):
        return general_params["dna_rna_lib_path"] if not self.lib_path else self.lib_path

class LipidLibrary(Library):
    def get_path(self):
        return general_params["lipid_lib_path"] if not self.lib_path else self.lib_path

class WaterLibrary(Library):
    def get_path(self):
        return general_params["water_lib_path"] if not self.lib_path else self.lib_path
