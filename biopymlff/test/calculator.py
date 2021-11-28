import unittest
import os
import time

from ase.io.proteindatabank import read_proteindatabank
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.mopac import MOPAC

import matplotlib.pyplot as plt

from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_gap import GEBF_GAP
from biopymlff.calculators.gebf_dft import GEBF_DFT

from biopymlff.train.ml import ML

class Calculator_Test(unittest.TestCase):
    def __init__(self, filename: str,  source_method: str, target_method: str, source: Calculator, target: Calculator):
        self.mol = read_proteindatabank(filename)
        self.source_method = source_method
        self.target_method = target_method
        self.source = source
        self.target = target

    def test_01_train(self):
        # Train the trainable calculators
        source_time = 0
        target_time = 0
        start_time = time.time()
        if isinstance(self.source, ML): self.source.train()
        end_time = time.time()
        source_time = end_time - start_time
        start_time = time.time()
        if isinstance(self.source, ML): self.target.train()
        end_time = time.time()
        target_time = end_time - start_time
        
        labels = []
        vals = []
        fig: plt.Figure = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])
        ax.bar(labels, vals)
        fig.savefig(os.getcwd() + "/benchmark/time/" + self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + ".train.png")

    def test_02_predict(self):
        control_pe = []
        experimental_pe = []
        
        control_time=[]
        experimental_time=[]

        samples = []
        self.mol.calc = MOPAC()
        dynamics = Langevin(atoms=self.mol, timestep=0.01, temperature_K=500, friction=0.01)
        collect_data = lambda: samples.append(atoms.copy())
        dynamics.attach(collect_data, interval=1)
        dynamics.run(steps=100)
        
        for sample in samples:
            time_start = time.time()
            source_pe = self.source.get_potential_energy()
            time_end = time.time()
            control_time.append(time_end - time_start)

            time_start = time.time()
            target_pe = self.target.get_potential_energy()
            time_end = time.time()
            experimental_time.append(time_end - time_start)
            
        plt.figure()
        plt.plot(control_pe, experimental_pe)
        plt.title(self.source_method 
            + " Potential Energy vs " 
            + self.target_method
            + " Potential Energy")
        plt.xlabel(self.source_method + " Potential Energy(J)")
        plt.ylabel(self.target_method + " Potential Energy(J)")
        plt.savefig(os.getcwd() 
            + "/benchmark/accuracy/" 
            +  self.source_method.replace(" ", "_").lower() 
            + "_pe_vs_" 
            + self.target_method.replace(" ", "_").lower() 
            + "_pe.png")
        
    
if __name__ == "__main__":
    unittest.main()