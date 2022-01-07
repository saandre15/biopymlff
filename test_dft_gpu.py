import os

from ase.io.proteindatabank import read_proteindatabank
from ase.io import write
from ase.md.langevin import Langevin

from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_pm6 import GEBF_PM6

from biopymlff.calculators.dft_gpu import DFT_GPU
from biopymlff.calculators.pm6_gpu import PM6_GPU


mol = read_proteindatabank(os.getcwd() + "/data/systems/proteins/4znn.pdb")
mol.calc = DFT_GPU(label="4znn_01_dft_gpu")

print(mol.get_potential_energy())
