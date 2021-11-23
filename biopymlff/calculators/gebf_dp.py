import os
import sys
import uuid
import random
import math

import numpy as np

from ase import Atoms
from ase.calculators.calculator import FileIOCalculator, ReadError, Parameters
from ase.units import kcal, mol, Debye
from ase.io import write
from ase.io.proteindatabank import read_proteindatabank

from ase.atoms import Atoms

from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin

from ase.calculators.gaussian import Gaussian
from ase.calculators.mopac import MOPAC

from quippy.potential import Potential
from quippy.descriptors import Descriptor

from deepmd.calculator import DP

import json

from biopymlff.calculators.gebf_ml import GEBF_ML
from biopymlff.helper import get_avaliable_gpu


class GEBF_DP(DP, GEBF_ML):
    """
    Generalized Energy Based Fragmentation utilitizing Deep Potential Machine Learning Technique

    Implemented Properties are `energy`, 'forces', and 'stress'

    Parameters
    ----------
    
    
    """

    datasplit = 0.5

    def __init__():
        
        super().__init__(restart, ignore_bad_restart_file, label, atoms, directory, pdb_id=pdb_id, ext_type="deep_pot")

    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):

        dataset_dir = self.data_dir + "/deep_pot_dataset"
        train_dir = dataset_dir + "/train"
        validate_dir = dataset_dir + "/validate"
        os.mkdir(dataset_dir)
        os.mkdir(train_dir)
        os.mkdir(validate_dir)

        set_index = 0
        # Convert this to numpy file
        for atom in traj:
            set_dir = (train_dir if random() < self.datasplit else validate_dir) + "/set_" + set_index
            # Figure this out
            box = atom.get_box()
            # Figure this out
            coord = atom.get_coord()
            forces = atom.get_forces()
            potential = atom.get_potential_energy()
            
            np.save(set_dir + "/box.npy", box)
            np.save(set_dir + "/coord.npy", coord)
            np.save(set_dir + "/forces.npy", forces)
            np.save(set_dir + "/potential.npy", potential)
            set_index+=1

        data = {}
        data['learning_rate'] = {
            'type': "exp",
            'start_lr': 0.001,
            'stop_lr': 3.51e-8,
            'decay_steps': 5000,
        }
        data['training'] = {}
        data['training']['training_data'] = {
            'systems': [train_dir],
            'batch_size': 'auto'
        }
        data['training']['validation_data'] = {
            "systems":		[validate_dir],
            "batch_size":	1,
            "numb_btch":	3
        }
        data['training'] = {
            "numb_step":	1000000,
            "seed":	 math.ceil(random() * 100000),
            "disp_file":	"lcurve.out",
            "disp_freq":	100,
            "save_freq":	1000
        }
        json_str = json.dumps(data)
        input_file = os.write("/tmp/" + uuid.uuid1() + ".json", json_str)
        
        if(len(get_avaliable_gpu()) != 0):
            os.environ['CUDA_VISIBLE_DEVICES'] = ""
            os.system("mpirun -l -launcher=fork -hosts=localhost -np 4 dp train --mpi-log=workers " + input_file)
        else:
            os.system("dp train " + input_file)
