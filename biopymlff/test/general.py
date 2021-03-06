
import unittest
import os
import time

from sklearn.metrics import mean_squared_error, r2_score

from ase.io.proteindatabank import read_proteindatabank
from ase.io import write
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.mopac import MOPAC
from ase.calculators.emt import EMT

from ase.md.langevin import Langevin

import matplotlib.pyplot as plt

import numpy as np

from ..train.ml import ML
from ..calculators.gebf_pm6 import GEBF_PM6
from ..calculators.gebf import GEBF
from ..cluster.slurm import run_parallel

class General_Test(unittest.TestCase):

    def setUp(self):
        if self.mol == None or self.source_method == None or self.target_method == None or self.source == None or self.target == None:
            raise NotImplementedError("Variables are not initalized")
        
    def init_variables(self, atom_filename: str, source_method: str, target_method: str, source: Calculator, target: Calculator):
        self.mol = read_proteindatabank(atom_filename)
        
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
        mol: Atoms = self.mol.copy()
        print("Atom Size " + str(len(mol)))
        mol.calc = GEBF_PM6(label="".join(mol.get_chemical_symbols()))
        dynamics = Langevin(atoms=mol, timestep=0.01, temperature_K=500, friction=1e-3, rng=np.random.default_rng(1000))
        collect_data = lambda: samples.append(mol.copy())
        dynamics.attach(collect_data, interval=1)
        dynamics.run(steps=100)

        write(os.getcwd() + "/data/traj/" + "".join(mol.get_chemical_symbols()), collect_data)
        
        # jobs = []
        samples = samples[0:len(samples) if len(samples) <= 50 else 50]

        for sample in samples:
            s: Atoms = sample
            time_start = time.time()
            source_mol = s.copy()
            source_mol.calc = self.source
            
            source_pe = source_mol.get_potential_energy()
            time_end = time.time()
            control_time.append(time_end - time_start)

            time_start = time.time()
            target_mol = s.copy()
            target_mol.calc = self.target
            target_pe = target_mol.get_potential_energy()
            time_end = time.time()
            experimental_time.append(time_end - time_start)
        
        # Plot out figure about potential energy
        rsme = mean_squared_error(control_pe, experimental_pe)
        r2 = r2_score(control_pe, experimental_pe)

        np.save(
            os.path.join(os.getcwd(), "benchmark", "accuracy", self.source_method + ".predict_pe.npy"), 
            np.array(control_pe)) # Save our control data file for future use
        np.save(
            os.path.join(os.getcwd(), "benchmark", "accuracy", self.target_method + "predict_pe.npy"),
            np.array(experimental_pe))# Save our experimental data file for future use
        # NOTE: Merge GPU and CPU results from different paritions

        fig = plt.figure()
        bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width*fig.dpi, bbox.height*fig.dpi
            
        plt.plot(control_pe, experimental_pe)
        plt.title(self.source_method 
            + " Potential Energy vs " 
            + self.target_method
            + " Potential Energy")
        plt.xlabel(self.source_method + " Potential Energy(J)")
        plt.ylabel(self.target_method + " Potential Energy(J)")
        plt.annotate("RMSE=" + rsme + "\nR^2=" + r2, ((width / 2) + 20, (height / 2) + 20))
        plt.savefig(os.getcwd() 
            + "/benchmark/accuracy/" 
            +  self.source_method.replace(" ", "_").lower() 
            + "_pe_vs_" 
            + self.target_method.replace(" ", "_").lower() 
            + "_pe.png")

        control_time_avg = sum(control_time) / len(control_time)
        experimiental_time_avg = sum(experimental_time) / len(control_time)

        # Plots out figure about average time it takes to perform a calculation
        fig = plt.figure()
        axes = fig.add_axes([0, 0, 1, 1])
        labels = [self.source_method, self.target_method]
        vals = [control_time_avg, experimiental_time_avg]
        plt.title("Time Average Between Methods(seconds))")
        axes.bar(labels, vals)
        plt.savefig(os.getcwd()
            + "/benchmark/time/"
            + self.source_method.replace(" ", "_").lower()
            + "_time_vs_"
            + self.target_method.replace(" ", "_").lower()
            + "_time.png")


