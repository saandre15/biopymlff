import os

from ase.io.proteindatabank import read_proteindatabank
from ase.io import write
from ase.md.langevin import Langevin

from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_pm6 import GEBF_PM6


mol = read_proteindatabank(os.getcwd() + "/data/systems/proteins/4znn.pdb")
mol.calc = GEBF_DFT(label="4znn_01")

print(mol.get_potential_energy())
