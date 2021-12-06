import os

from ase.io.proteindatabank import read_proteindatabank
from ase.md.langevin import Langevin

from biopymlff.calculators.gebf_dft import GEBF_DFT


mol = read_proteindatabank(os.getcwd() + "/data/systems/proteins/4znn.pdb")
mol.calc = GEBF_DFT(label="4znn_01")

print(mol.get_potential_energy())