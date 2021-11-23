import os

from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.md.langevin import Langevin

from deepmd.calculator import DP

from ase.io import write
from ase.io.proteindatabank import read_proteindatabank

import numpy as np

from biopymlff.train.library import LibraryType

class ML():
    """
    Machine Learning class to extend Calculators with machine learning capabilities

    """
    
    def __init__(self, descriptors: np.ndarray, lib_type: LibraryType, lib_path: str=None):
        self.descriptors = descriptors
        self.models = {}
        if lib_predef:
            if lib_path_or_type == "protein":
                pass
            elif lib_path_or_type == "dna+rna":
                pass
            elif lib_path_or_type == "lipids":
                pass
            elif lib_path_or_type == "water":
                pass
        self.lib_path = lib_path_or_type
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

    def add_to_lib(self, atoms: Atoms):
        self.lib_cache.append(atoms)
        if not os.path.exists(self.lib_path): os.mkdir(self.lib_path)
        
        write(self.lib_path + "/" )

    def load_lib(self) -> list:
        files = [file for file in os.listdir(self.lib_path) if os.path.isfile(self.lib_path + "/" + file)]
        atoms = list(np.zeros(shape=(len(files), 1)))
        counter = 0
        for file in files:
            atom = read_proteindatabank(file)
            atoms[counter] = atom
            counter+=1
        return atoms
            
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