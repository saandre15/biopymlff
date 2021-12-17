import unittest
import os
import time
import math

from ase.io.proteindatabank import read_proteindatabank
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.md.langevin import Langevin

# from pytraj.analysis import dihedral_analysis

import matplotlib.pyplot as plt

from sklearn.metrics import mean_squared_error, r2_score

import numpy as np

from ..test.biomol import Biomol_Test
from ..util.getenv import getenv
from ..util.convert import ase_atoms_to_pytraj_atoms
from ..train.library import ProteinLibrary
from ..train.ml import ML


class Protein_Test(Biomol_Test):
    
    def __init__(self, ):
        super.__init__()
        

    def test_01_absolute_deviation(self): # Run this on GEBF-DP, ANI-1x and ff14SB and merge all the graph to replicate Fig. 7 on GEBF_ML papers
        source_traj = []
        target_traj = []
        conformer_indexes = [] # TODO: Figure out how to randomly sample base off a seed

        fig = plt.figure()
        plt.title("Absolute Deviation by Conformers")
        ax = fig.add_axes([0, 0, 1, 1])
        
        deviations_by_conformers = []

        for index in conformer_indexes:
            source = source_traj[index]
            target = target_traj[index]
            
            source_energy = source.get_potential_energy()
            target_energy = target.get_potential_energy()

            deviation = abs(target_energy - source_energy)
            deviations_by_conformers.append(deviation)
        
        X = np.arange(len(deviations_by_conformers))
        ax.bar(X + 0.00, data, color='b', width=0.25)
        plt.savefig(os.path.join(
            os.getcwd(),
            "benchmark",
            "accuracy",
            self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + "absolute_deviation.png"
        ))
            

    def test_02_dihedral_anaylsis(self):
        
        source_traj = []
        target_traj = []

        analysis = {
            '\psi': [dihedral_analysis.calc_psi(traj=source_traj, range360=True), dihedral_analysis.calc_psi(traj=target_traj, range360=True)],
            '\omega': [dihedral_analysis.calc_omega(traj=source_traj, range360=True), dihedral_analysis.calc_omega(traj=target_traj, range360=True)],
            '\phi': [dihedral_analysis.calc_phi(traj=source_traj, range360=True), dihedral_analysis.calc_phi(traj=self.target_traj, range360=True)]
        }

        for key in analysis:
            x_label = key
            y_label = "Distribution (10^-3)"
            source=analysis[key][0]
            target=analysis[key][1]
            fig = plt.figure()
            plt.title(key)
            plt.plot(source, label=self.source_method)
            plt.plot(target, label=self.target_method)
            os.path.join(os.getcwd(), 
                "benchmark", 
                "accuracy", 
                self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + "dihedral_analysis.png")

    

    def test_03_end_to_end(self):
        G = AtomGraph(self.mol)
        start_atom, end_atom = G.furthest_atoms()
        start_index = G.find(start_atom)
        end_index = G.find(end_atom)

        source_traj = []
        target_traj = []
        
        dist = dihedral_analysis.distance(source_traj, [[start_index], [end_index]])

        fig = plt.figure()
        plt.plot(dist)
        plt.xlabel("t(ns)")
        plt.ylabel("distance(A)")
        
        plt.save
        
    
    def test_04_rsme_energy_and_force(self):
        # TODO: Figure out how to train mulitple iterations with different dataset seize

        if not isinstance(self.source, ML): return; # Skips this test if the calculator is not trainable


        train_size = [1000, 2000, 3000, 4000, 5000] # Dataset train size

        for size in train_size:
            lib = ProteinLibrary()
            datasets = lib.extract_atoms(size) # Randomly select atoms that output 
            
            self.source.reset()
            for dataset in datasets:
                self.source.train()
            
        

        # As ASE Atoms
        source_traj = []
        target_traj = []
        
        source_count = len(source_traj)
        target_count = len(target_traj)
        
        source_pe = [val.get_potential_energy() for val in source_traj]
        target_pe = [val.get_potential_energy() for val in target_traj]
        
        source_forces = [np.ndarray.item(val.get_forces()) for val in source_traj]
        target_forces = [np.ndarray.item(val.get_forces()) for val in target_traj]
        
        rsme_E = mean_squared_error(target_pe, source_pe)
        rsme_forces = mean_squared_error() # convert forces to scalars

        fig = plt.figure()
        plt.title("RSME")
        plt.plot(source, label=self.source_method)
        plt.plot(target, label=self.target_method)
        os.path.join(os.getcwd(), 
            "benchmark", 
            "accuracy", 
            self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + ".png")
        

        