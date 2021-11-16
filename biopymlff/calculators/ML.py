import os

from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.md.langevin import Langevin

class ML(Calculator):

    _deprecated = object()
    
    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.', pdb_id=None,
                 **kwargs):

        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, kwargs=kwargs)

        self.models = {}
        self.pdb_id = pdb_id
        self.data_dir=os.getcwd() + "/data/" + self.pdb_id


    # Validates if the models have bee
    def add_model(self, type: str, model_path: str):
        self.models[type] = model_path
    
    def get_model(type: str):
        return self.models[type]

    def calculate(self, atoms: Atoms):
        for model in self.models:
            if not os.path.exists(model): raise IOError("Model not found.")
    
    def train(self, atoms: Atoms):
        raise NotImplementedError("train has not been implemented.")

    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):
        raise NotImplementedError("train_model has not been implemented.")

    def generate_subsets(self, atoms: Atoms, calc: Calculator):      
        temp = 500
        friction = 0.1
        
        atoms.calc = calc

        traj_db=self.run_md(atoms, 0.05, temp, friction)
        
        return traj_db
        
    def run_md(self, atoms: Atoms, timestep: float, temp: float, friction: float):
        print("Running MD Simulation")
        traj_db = []
        dynamics = Langevin(atoms=atoms, timestep=timestep, temperature_K=temp, friction=friction)
        collect_data = lambda: traj_db.append(atoms.copy())
        dynamics.attach(collect_data, interval=10)
        dynamics.run(steps=100)
        
        return traj_db