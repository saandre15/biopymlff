import random
import math
import os

from ase.atoms import Atoms
from ase.atom import Atom
from ase.geometry.analysis import Analysis
from ase.io.proteindatabank import read_proteindatabank

from ase.io import write

from biopymlff.data.AtomGraph import AtomGraph

symbolA = 'C'
symbolB = 'C'

mol = read_proteindatabank(os.getcwd() + "/data/systems/4znn+h.pdb")
write(filename=os.getcwd() + "/data/systems/4znn+h.mol2", images=mol)