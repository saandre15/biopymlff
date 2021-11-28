import os

from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.md.langevin import Langevin

from ase.io import write
from ase.io.proteindatabank import read_proteindatabank

import numpy as np

from ..descriptors.descriptor import Descriptor
from ..train.library import Library

class ML():
    """
    Machine Learning class to extend Calculators with machine learning capabilities

    """
    
    def __init__(self, descriptors: Descriptor, library: Library):
        self.descriptors = descriptors
        self.models = {}
        if lib_path == None:
            lib_path = lib_type.get_path()
        self.lib_path = lib_path
        self.lib_cache = self.load_lib()


    def add_model(self, type: str, model_path: str):
        """
        Registers ML model files

        Notes
        -----
        Used to check the model file if the model file exists and if not create one or train an existing model.

        Parameters
        ----------
        type: str
            Model ID
        model_path
            Path to to the model file
        """
        self.models[type] = model_path
    
    def get_model(self, type: str) -> str:
        """
        Returns the stored model.

        Parameters
        ---------- 
        type: str
            Model ID
        
        Returns
        -------
        Returns the model file path
        """
        return self.models[type]


    def train(self, atoms: Atoms):
        """
        Trains all of the saved models using different libraries or self implemented algorithms to train.
        
        Parameters
        ----------
        atoms: Atoms 
            ASE Atoms
        """
        raise NotImplementedError("train has not been implemented.")

    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):
        """
        Train each individual models 

        Parameters
        ----------
        model_file: str
            Path to to the model file
        atypes: list
            Atom Symbol List within the Trajectory
        traj: list
            List of ASE Atoms in a timeframe list format
        type: str
            Model ID
        """
        raise NotImplementedError("train_model has not been implemented.")

    def generate_subsets(self, atoms: Atoms, calc: Calculator) -> list:
        """
        Creates subset of ASE Atoms as training or validation dataset for the ML program or algorithms

        Parameters
        ----------
        atoms: Atoms
            ASE Atoms
        calc: Calculator
            ASE Calculators used to extract out the potential energy for MD simulation
    
        Returns
        -------
        List of ASE Atoms that have been permutated by MD simulation.
        """      
        temp = 500
        friction = 0.1
        
        atoms.calc = calc

        traj_db=self.run_md(atoms, 0.05, temp, friction)
        
        return traj_db
        
    def run_md(self, atoms: Atoms, timestep: float, temp: float, friction: float) -> list:
        """
        Runs MD simulation to extract out snapshots of each configuration to be processed.

        Parameters
        ----------
        atoms: Atoms
            Atoms 
        timestep: float
            Timestep for numerical analysis calculations
        temp: float
            Tempature at which the MD simulation will be runned at.
        friction: float
            
        Returns
        -------
        List of atoms trajectory         
        """
        print("Performing MD simulation to extract out dataset.")
        traj_db = []
        dynamics = Langevin(atoms=atoms, timestep=timestep, temperature_K=temp, friction=friction)
        collect_data = lambda: traj_db.append(atoms.copy())
        dynamics.attach(collect_data, interval=10)
        dynamics.run(steps=100)
        print("MD simulation completed.")
        
        return traj_db