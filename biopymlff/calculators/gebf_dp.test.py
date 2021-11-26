import unittest
import os

from ase.io.proteindatabank import read_proteindatabank

from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_gap import GEBF_GAP


class GEBF_GAP_Test(unittest.TestCase):

    def __init__():
        # super(self, unittest.TestCase).__init__()
        self.lib = ProteinLibrary()
        self.mol = read_proteindatabank(os.getcwd() + "/data/4znn.pdb")
        self.calc = GEBF_GAP(directory=os.getcwd() + "/data/4znn/GEBF_GAP/", library=self.lib, atoms=mol)

    def test_01_train(self):
        print("Started training")
        self.calc
        print("Ended training")

    def test_02_predict(self):
        self.calc = e        
        

if __name__ == "__main__":
    unittest.main()