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
from biopymlff.calculators.GEBF_DFT import GEBF_DFT
from biopymlff.calculators.guassian import Gaussian


class GEBF_ML(GEBF_PM6, ML):

    implemented_properties = ['energy', 'energies', 'forces', 'stresses']

    _deprecated=Calculator._deprecated

    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.', pdb_id=None, ext_type=None,
                 **kwargs):

        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, directory=directory, pdb_id=pdb_id)

        self.ext_type = ext_type
        self.add_model("dft", self.dft_model_file)
        self.add_model("pm6", self.pm6_model_file)

        self.sigmas = []
        self.epsilon_bayes = 0

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
        subsystems=self.subsystems()

        config=toml.load(os.getcwd() + "/env.toml")

        g_config=config["gaussian"]

        for subsys in subsystems:
            atoms: Atoms = subsys

            pm6_traj = self \
                .generate_subsets(atoms, MOPAC(
                    method="PM6"
                ))
            
            for index in range(0, len(pm6_traj)):
                pm6_snapshot = pm6_traj[index]
                updated=False
                sigma_i = self.calculate_sigma()
                steps_after_qm = 0
                m = 10
                self.sigmas.append(sigma_i)
                sigma_max = max(self.sigmas)
                # Decides training
                if sigma_max > 2 * self.epsilon_bayes:
                    if sigma_i > self.epsilon_bayes:
                        self.update_model(pm6_snapshot)
                        updated = True
                        steps_after_qm+=1
                        m = max(10, 0.18 / self.calculate_max_force_error(pm6_snapshot))
                    else:
                        if self.in_lib(pm6_snapshot):
                            updated = False
                        else:
                            self.update_model(pm6_snapshot)
                            updated = True
                            steps_after_qm+=1
                            m = max(10, 0.18 / self.calculate_max_force_error(pm6_snapshot))
                else:
                    if steps_after_qm < m:
                        if sigma_max > self.epsilon_bayes:
                            if sigma_i > self.epsilon_bayes:
                                self.update_model(pm6_snapshot)
                                updated = True
                                steps_after_qm+=1
                                m = max(10, 0.18 / self.calculate_max_force_error(pm6_snapshot))
                            else:
                                if self.in_lib(pm6_snapshot):
                                    updated = False
                                else:
                                    self.update_model(pm6_snapshot)
                                    updated = True
                                    steps_after_qm+=1
                                    m = max(10, 0.18 / self.calculate_max_force_error(pm6_snapshot))
                    else:               
                        updated = False

                # Decides epsilon_bayes
                if updated:
                    sigmas_after_qm = self.sigmas[len(self.sigmas) - 1 - steps_after_qm:len(self.sigmas) - 1]
                    sigma_max = max(sigmas_after_qm)
                    if(steps_after_qm > 9):
                        if self.get_standard_deviation(sigmas_after_qm) < 0.2:
                            self.epsilon_bayes = (1/10) * sum(sigmas_after_qm)            
                
    def update_model(self, atoms: Atoms):
        symbols = list(set(atoms.get_chemical_symbols()))
        self.train_model(self.models("dft"), symbols, [atoms], "dft")
        self.train_model(self.models("pm6"), symbols, [atoms], "pm6")

    def calculate_sigma(self):
        pass

    def calculate_max_force_error(self, atoms: Atoms):
        pass

    def in_lib(self, atoms: Atoms) -> list:
        pass

    def get_gaussian(self, method=None, basis=None):
        general_params = getenv()['general']
        gaussian_params = getenv()['gaussian']
        return Gaussian(
            label=self.label,
            directory=general_params['scratch_dir'] + "/gaussian/" + self.label + "_" + uuid.uuid1().hex,
            mem=gaussian_params['mem'] \
                if gaussian_params['mem'] != "auto" \
                else str((os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024.**3)) - 5) + "GB",
            chk=gaussian_params['checkpoint_file'] if gaussian_params['checkpoint_file'] else self.label + ".chk",
            save=None,
            method=method,
            basis=basis,
            scf='(noincfock,novaracc,fermi,maxcycle=3000,ndamp=64,xqc)',
            int='acc2e=12'
        )
       
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

        