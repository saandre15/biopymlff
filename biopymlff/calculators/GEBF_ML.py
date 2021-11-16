import os
import sys
import uuid
import random
import toml

import numpy as np

from ase import Atoms, Atom
from ase.calculators.calculator import Calculator, ReadError, Parameters
from ase.units import kcal, mol, Debye
from ase.io import write
from ase.io.proteindatabank import read_proteindatabank
from ase.io.xyz import read_xyz

from ase.atoms import Atoms

from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin

from ase.calculators.mopac import MOPAC

from quippy.potential import Potential

from biopymlff.calculators.ML import ML
from biopymlff.calculators.GEBF_PM6 import GEBF_PM6
from biopymlff.calculators.guassian import Gaussian


class GEBF_ML(GEBF, ML):

    implemented_properties = ['energy', 'energies', 'forces', 'stresses']

    _deprecated=Calculator._deprecated

    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.', pdb_id=None, ext_type=None,
                 **kwargs):

        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, directory=directory, pdb_id=pdb_id)

        self.ext_type = ext_type
        self.add_model("dft", self.dft_model_file)
        self.add_model("pm6", self.pm6_model_file)

    def get_subfrag_dir(self): return self.data_dir + "/" + self.pdb_id + "_subsys"

    def calculate(self, atoms: Atoms):
        gap_dft_pot=Potential(param_filename=self.get_model("dft"))
        gap_pm6_pot=Potential(param_filename=self.get_model("pm6"))
        
        gebf_pm6_pot = GEBF_PM6()

        atoms.calc = gebf_pm6_pot
        
        pm6_base_energy = atoms.get_potential_energy()
        
        atoms.calc = gap_dft_pot
        dft_energy=atoms.get_potential_energy()
        atoms.calc = gap_pm6_pot
        pm6_energy=atoms.get_potential_energy()
        
        energy = pm6_base_energy + (dft_energy - pm6_energy)

        self.results = {
            'energy': energy,
            'forces': self.calculate_numerical_forces(atoms),
            'stresses': self.calculate_numerical_stress(atoms)
        }

    def train(self, atoms: Atoms):
        # Runs some kind of structal optimization if neccesiary [MOPAC]
        subsystems=self.subfrag(atoms)
            
        all_dft_traj = []
        all_pm6_traj = []
        atom_types = []

        config=toml.load(os.getcwd() + "/env.toml")

        g_config=config["gaussian"]

        for subsys in subsystems:
            atoms = subsys

            dft_traj = self \
                .generate_subsets(atoms, Gaussian(
                    mem=g_config["memory"],
                    chk=g_config["checkpoint_file"],
                    save=None,
                    method=g_config["method"],
                    basis=g_config["basis"],
                    scf="qc"
                ))

            all_dft_traj+=(dft_traj)

            pm6_traj = self \
                .generate_subsets(atoms, MOPAC(
                    method="PM6"
                ))

            all_pm6_traj+=(pm6_traj)

            atom_types+=atoms.get_chemical_symbols()

        atom_types = list(set(atom_types))            
            
        # Figure out some sort of way to merge
        
        self.train_model(self.dft_model_file, atom_types, all_dft_traj)
        self.train_model(self.pm6_model_file, atom_types, all_pm6_traj)


    # See FIG S7
    def descriminate(self, atoms: Atoms) -> bool:
        pass

    def parse_fragment(self, line: str, system: Atoms) -> Atoms:
        fields = line.split()
        index = fields[0]
        unknown_a = fields[1]
        subfrags_as_str = fields[2]
        unknown_b = fields[3]
        
        subfrags_as_array = subfrags_as_str \
            .replace(")", "") \
            .replace("(", "") \
            .split(",") \
        
        atoms_symbol = system.get_chemical_symbols()
        atoms_pos = system.get_positions()
        atoms = []

        for subfrag in subfrags_as_array:
            try:
                start = int(subfrag.split("-")[0]) - 1
                end = int(subfrag.split("-")[1]) if len(subfrag.split("-")) > 1 else start + 1
                _range = end - start
                # print("start " + str(start))
                # print("end " + str(end))
                # print("range " + str(_range))
                for index in range(_range):
                    atom = Atom(symbol=atoms_symbol[index + start], position=atoms_pos[index + start])
                    atoms.append(atom)
            except IndexError:
                pass
                
        mol = Atoms(atoms)
        return mol

        