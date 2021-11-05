import math
import numpy as np
import scipy as sp

from ase.atoms import Atoms

import tensorflow as tf
from tensorflow.keras.models import load_model, Model

from mendeleev import element

from biopymlff.calculators.ML import ML


class GEBF_HDNN4G(ML):

    implemented_properties = ['energy', 'energies', 'force', 'stresses']

    def __init__(self):

        self.atomic_electronegativity_model_file=self.data_dir + "/electronegativity.hdnn4g.h5"
        self.atomic_shortrange_energy_model_file=self.data_dir + "/short_range_energy.hdnn4g.h5"
        self.add_model(self.atomic_electronegativity_model_file)
        self.add_model(self.atomic_shortrange_energy_model_file)

    def calculate(self, atoms: Atoms):
        super().calculate(atoms)
        
        electronegativity_model: Model = load_model(self.atomic_electronegativity_model_file)
        shortrange_energy_model: Model = load_model(self.atomic_shortrange_energy_model_file)

        colomb_energy = 0
        short_energy = 0

        symbols = atoms.symbols()
        size = len(symbols)
        coords = atoms.get_positions()
        charges = atoms.get_initial_charges()
        hardness = []
        # Testing with hardness pre charge modification
        for symbol in symbols:
            el = element(symbol)
            el.hardness()
        charge_den_of_width =[]

        short_interatomic_potential = []
        elec_interatomic_potential = []
        interatomic_potential = []
        
        sym_func = []
        
        for desc in sym_func:
            val = electronegativity_model.predict(desc)
            electronegativities.append(val[0])
        
        charges = self.equilibrate_charges(coords, charges, electronegativities, hardness, charge_den_of_width)

        electronegativities = []

        # Testing with hardness post charge modification
        # for symbol in symbols:
        #     el = element(symbol)
        #     el.hardness()
        
        for index in size:
            desc_1 = charges[index]
            desc_2 = sym_func[index]
            short_interatomic_potential.append(neural_network.predict())
            elec_interatomic_potential.append(
                self.get_interatomic_elec_energy(
                    coords[index], 
                    coords, 
                    cur_charge, 
                    charges, 
                    cur_hardness, 
                    hardnesses, 
                    charge_den_of_width, 
                    electronegativity
                )
            )
            interatomic_potential.append(short_interatomic_potential[index] + elec_interatomic_potential[index])


        self.results = {
            'energy': sum(short_interatomic_potential) + sum(elec_interatomic_potential),
            'energies': interatomic_potential,
            'forces': self.calculate_numerical_forces(atoms),
            'stresses': self.calculate_numerical_stress(atoms)
        }

    
    def get_interatomic_elec_energy(self, 
        cur_coord: float, 
        coords: list, 
        cur_charge: float, 
        charges: list, 
        cur_hardness: float,
        hardnesses: list, 
        charge_den_of_width: float, 
        electronegativity: float
    ):
        return self.get_interatomic_ref_elec_energy(
            cur_coord, coords, cur_charge, charges, cur_hardness, hardnesses, electronegativity
        ) + electronegativity * cur_charge + 0.5 * cur_hardness * cur_charge * cur_charge

    def get_interatomic_ref_elec_energy(
        cur_coord: float, 
        coords: list, 
        cur_charge: float, 
        charges: list, 
        cur_hardness: float,
        hardnesses: list, 
        charge_den_of_width: float, 
        electronegativity: float
    ):
        size = len(coords)
        total = 0
        for index in size:
            next_coord = coords[index]
            next_charge = charges[index]
            next_hardness = hardnesses[index]
            radius = self.radius(cur_coords, next_coords)
            total+=(math.erf(radius / (math.sqrt(2) * math.sqrt(cur_hardness * cur_hardness + next_hardness + next_hardness))) * cur_charge * next_charge) / radius
            
        total+=(charge * charge) / (2 * charge_den_of_width * math.sqrt(math.pi))
        return total

    # Looks for the minimum potential energy contributed from the charge while also keeping the same system total charge the same.
    # Based Off CENT
    def equilibrate_charges(self, coord: list, charges: list, electronegativities: list, hardness: list, charge_den_of_width: list):
        total=0
        for charge in charges:
            total+=charge
        # Makes sure the charges is always equal to the original system charges
        constraints = ({ 'type': 'eq', 'fun': lambda func: sum(charges) - total})
        # Does this create a seperate charge list to perform calculations or does it edit the orignal?
        new_charges = sp.optimize.minimize(
            lambda x: self.get_elec_energy(x[0], x[1], x[2], x[3], x[4]), [[coord, charges, electronegativities, hardness, charge_den_of_width]], 
            method='COBYLA', 
            constraints=constraints)
        return new_charges

    def radius(a: list, b: list):
        x_diff = a[0] - b[0]
        y_diff = a[1] - b[1]
        z_diff = a[2] - b[2]
        return math.sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff) 
     
    
    def train(self, atoms: Atoms):
        
        dataset_dir = self.data_dir + "/pm6_dataset"

        

        self.generate_subsets(atoms, )
        

    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):
        if type == "default":
            raise ValueError("Make sure the model")
        elif type == "electronegativity":
            pass
        elif type == ""

        


    
