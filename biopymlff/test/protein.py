import unittest
import os
import time

from ase.io.proteindatabank import read_proteindatabank
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator

import matplotlib.pyplot as plt

from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_gap import GEBF_GAP
from biopymlff.calculators.gebf_dft import GEBF_DFT

from biopymlff.train.ml import ML

class Protein_Test(unittest.TestCase):

    def __init__(self):
        pass
    
    def test_01_absolute_deviation(self):
        pass

    def test_02_dihedral_anaylsis(self):
        pass

    def test_03_end_to_end(self):
        pass
    
    def test_04_mae_force(self):
        pass

    def test_05_rmsd_from_init(self):
        pass
    
    def test_06_rsme_energy_and_force(self):
        pass