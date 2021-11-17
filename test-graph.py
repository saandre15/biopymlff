import os
import random
import math

from ase.atoms import Atoms
from ase.atom import Atom
from ase.geometry.analysis import Analysis
from ase.io.proteindatabank import read_proteindatabank

from biopymlff.data.AtomGraph import AtomGraph

symbolA = 'C'
symbolB = 'C'

mol = read_proteindatabank(os.getcwd() + "/data/systems/4znn+h.pdb")
graph = AtomGraph(mol)
print("SIZE " + str(graph.size()))
graph.show()
# fragments = graph.fragments_by_bond_as_indexes(symbolA, symbolB)
# print(fragments)
atoms = graph.fragment_by_bond_as_atoms_list(symbolA, symbolB)
print(atoms)

counter = 0
for atom in atoms:
    g = AtomGraph(atom)
    # if g.size() == 1: continue
    # counter+=len(atom.get_positions())
    counter+=g.size()
    print("Counter " + str(counter))
    g.show()