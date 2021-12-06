import unittest
import os
import time

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


class Protein_Test(Biomol_Test):

    def test_01_absolute_deviation(self):
        # TODO: Figureout conformers?
        pass

    def test_02_dihedral_anaylsis(self):
        pass
        
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
            plt.savefig(os.getcwd() 
                + "/benchmark/accuracy/" 
                + self.source_method.replace(" ", "_").lower() 
                + "_vs_" 
                + self.target_method.replace(" ", "_").lower())

    

    def test_03_end_to_end(self):
        pass
        G = AtomGraph(self.mol)
        start_atom, end_atom = G.furthest_atoms()
        start_index = G.find(start_atom)
        end_index = G.find(end_atom)

        source_traj = []
        target_traj = []

        index_one = -1
        index_two = -1
        index_three = -1
        

        dist = dihedral_analysis.distance(source_traj, [[start_index], [end_index]])

        fig = plt.figure()
        plt.plot(dist)
        plt.xlabel("t(ns)")
        plt.ylabel("distance(A)")


    # def test_04_mae_force(self):
    #     pass

    # def test_05_rmsd_from_init(self):
    #     pass
    
    def test_06_rsme_energy_and_force(self):
        pass
        # TODO: Figure out how to train mulitple iterations with different dataset seize

        # mol = self.mol.copy()

        # # As ASE Atoms
        # source_traj = []
        # target_traj = []
        
        # source_count = len(source_traj)
        # target_count = len(target_traj)
        
        # source_pe = [val.get_potential_energy() for val in source_traj]
        # target_pe = [val.get_potential_energy() for val in target_traj]
        
        # source_forces = [np.ndarray.item(val.get_forces()) for val in source_traj]
        # target_forces = [np.ndarray.item(val.get_forces()) for val in target_traj]
        
        # rsme_E = mean_squared_error(target_pe, source_pe)
        