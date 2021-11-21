import os

from ase.io.proteindatabank import read_proteindatabank

from biopymlff.calculators.GEBF import GEBF

mol = read_proteindatabank(os.getcwd() + "/data/systems/4znn.pdb")

print(mol)
gebf = GEBF(label="4znn", atoms=mol)
mol.calc = gebf
mol.get_potential_energy()