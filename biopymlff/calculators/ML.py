import os

from ase.atoms import Atoms
from ase.calculators.calculator import Calculator

class ML(Calculator):

    _deprecated = object()
    
    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.',
                 **kwargs):

        self.models = []
        self.data_dir=os.getcwd() + "/data/" + self.pdb_id
        
        # Checks if training is required

    # Validates if the models have bee
    def add_model(self, model: str):
        self.models.append(model)

    def calculate(self, atoms: Atoms):
        for model in self.models:
            if not os.path.exists(model): raise IOError("Model not found.")
    
    def train(self, atoms: Atoms):
        raise NotImplementedError("train_model has not been implemented.")

    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):
        raise NotImplementedError("train_model has not been implemented.")


    def generate_subsets(self, atoms: Atoms, calc: Calculator):      
        temp = 500
        friction = 0.1

        traj_db=self.run_md(atoms, desc_vector, temp, friction)
        
        return traj_db
        
    def run_md(self):
        traj_db = []
        
        dynamics = Langevin(atoms, timestep)
        collect_data = lambda: traj_db.append(atoms.copy())
        dynamics.attach(collect_data, interval=10)
        dynamics.run(steps=500)
        
        return traj_db