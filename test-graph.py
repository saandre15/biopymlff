import os
import random
import math

from ase.atoms import Atoms
from ase.atom import Atom
from ase.geometry.analysis import Analysis
from ase.io.proteindatabank import read_proteindatabank

from biopymlff.util.AtomGraph import AtomGraph

symbolA = 'C'
symbolB = 'C'

mol = read_proteindatabank(os.getcwd() + "/data/systems/4znn.pdb")
graph = AtomGraph(mol)
graph.fragments_by_bond(symbol, symbolB)
