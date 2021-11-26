import unittest
import os
import time

from ase.io.proteindatabank import read_proteindatabank
from ase.atoms import Atoms

import matplotlib.pyplot as plt


from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_gap import GEBF_GAP
from biopymlff.calculators.gebf_dft import GEBF_DFT


class GEBF_GAP_Test(unittest.TestCase):

    def __init__():
        # super(self, unittest.TestCase).__init__()
        self.lib = ProteinLibrary()
        self.mol = read_proteindatabank(os.getcwd() + "/data/4znn.pdb")
        self.calc = GEBF_GAP(directory=os.getcwd() + "/data/4znn/GEBF_GAP/", library=self.lib, atoms=mol)

    def test_01_train(self):
        start_time = time.time()
        self.calc.train(self.mol)
        end_time = time.time()
        os.write(os.getcwd() + "/benchmark/time/gebf_gap.train.txt", str(end_time - start_time) + " seconds")

    def test_02_predict(self):
        # DFT Energy
        control = []
        # GEBF_GAP Energy
        experimental = []
        gebf_gap: GEBF_GAP = self.calc
        samples = gebf_gap.run_md(self.mol, 0.1, 500, 0.1)
        
        for sample in samples:
            s: Atoms = sample
            s.calc = gebf_gap
            gebf_gap_energy = s.get_potential_energy()
            s.calc = GEBF_DFT()
            dft_energy = s.get_potential_energy()
            control.append(dft_energy)
            experimental.append(gebf_gap_energy)
        
        plt.plot(control, experimental)
        plt.title("DFT PE vs GEBF_GAP PE")
        plt.xlabel("DFT Potential Energy(Joule)")
        plt.ylabel("GEBF_GAP Potential Energy(Joule)")
        plt.savefig(os.getcwd() + "/benchmark/accuracy/gebf_dft_pe_vs_gebf_gap_pe.png")
        
    
if __name__ == "__main__":
    unittest.main()