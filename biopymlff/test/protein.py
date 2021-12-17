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
        

    def test_01_absolute_deviation(self): # Run this on GEBF-DP, ANI-1x and ff14SB and merge all the graph to replicate Fig. 7 on GEBF_ML papers
        source_traj = self.get_source_frames()
        target_traj = self.get_target_frames()
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
            plt.savefig(os.path.join(os.getcwd(), 
                "benchmark", 
                "accuracy", 
                self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + "dihedral_analysis.png"))

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
        
        plt.savefig(os.path.join(os.getcwd(), "benchmark", "accuracy", 
            self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + "end_to_end.png"))
        
    
    def test_04_rsme_energy_and_force(self):

        if not isinstance(self.source, ML): return; # Skips this test if the calculator is not trainable

        train_size = [1000, 2000, 3000, 4000, 5000] # Dataset train size
        rsme_E = []
        rsme_F = []

        for size in train_size:
            lib = ProteinLibrary()
            datasets = lib.extract_atoms(size) # Randomly select atoms that output subsystem of size N
            
            self.source.reset()
            for dataset in datasets:
                self.source.train(dataset)
            
            traj = []
            pes = []
            forces = []
            
            for atom in traj:
                pes.append(atom.get_potential_energy())
                for forces in atom.get_forces():
                    force_scalars = forces
                    forces.append(atom.get_forces()) # NOTE: Convert to scalars

            # Continue using RMSE E and F

        fig = plt.figure() # Energy RSME vs Subsystem Dataset Size Bar Graphs
        plt.title("Energy RSME vs Subsystem Dataset Size")
        ax = fig.add_axes([0, 0, 1, 1])
        X = np.arange(len(train_size))
        ax.bar(X + 0.00, train_size, rsme_E, color='b', width=0.25)
        plt.savefig(os.path.join(os.getcwd(), 
            "benchmark", 
            "accuracy", 
            self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + ".energy_rmse_vs_train_dataset_size.png"))

        fig = plt.figure() # Force RSME vs Subsystem Dataset Size Bar Graphs
        plt.title("Force RSME vs Subsystem Dataset Size")
        ax = fig.add_axes([0, 0, 1, 1])
        X = np.arange(len(train_size))
        ax.bar(X + 0.00, train_size, rsme_F, color='b', width=0.25)
        plt.savefig(os.path.join(os.getcwd(), 
            "benchmark",
            "accuracy",
            self.source_method.replace(" ", "_").lower() + "_vs_" + self.target_method.replace(" ", "_").lower() + "forces_rmse_vs_train_dataset_size.png"))

        
        

        