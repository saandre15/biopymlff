import os
import sys
import uuid
import random

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

from ase.calculators.gaussian import Gaussian
from ase.calculators.mopac import MOPAC

from quippy.potential import Potential
from quippy.descriptors import Descriptor

class GEBF_ML(Calculator):

    implemented_properties = ['energy', 'energies', 'forces', 'stresses']

    _deprecated=object()

    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.',
                 **kwargs):

        self.pdb_id = kwargs.get("pdb_id")
        self.ext_type = kwargs.get("ext_type")
        
        self.data_dir=os.getcwd() + "/data/" + self.pdb_id
        self.subfrag_dir=self.data_dir + "/" + self.pdb_id + "_subsys"
        self.dft_model_file=self.data_dir + "/dft_model.{self.ext_type}.xml"
        self.pm6_model_file=self.data_dir + "/pm6_model.{self.ext_type}.xml"

        
        # Verifies if the models has been created
        if (not os.path.exists(self.data_dir) \
            or \
            (not os.path.exists(self.dft_model_file)) and (not os.path.exists(self.pm6_model_file))):
            self.train(atoms)
        
        Calculator.__init__(self, restart, ignore_bad_restart_file, label, atoms, kwargs)

    def calculate(self, atoms: Atoms):
        gap_dft_pot=Potential(param_filename=self.dft_model_file)
        gap_pm6_pot=Potential(param_filename=self.pm6_model_file)
        
        mopac_pot = Potential(
            calculator=MOPAC(
                method="PM6"
            )
        )

        atoms.calc = mopac_pot
        
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
        temp = 300
        friction = 1
        # Runs some kind of structal optimization [MOPAC]
        subfrag_dir = self.subfrag(atoms)
        filename_woext=""
        numfiles = len([f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f)) and f[0] != '.'])
        
        frag_count = numfiles

        all_dft_traj = []
        all_pm6_traj = []
        atom_types = []

        for index in range(1, frag_count):
            atoms: Atoms = read_xyz(subfrag_dir + filename_woext + "_" + index + ".xyz")
            
            dft_traj = self \
                .generate_subsets(atoms, Potential(
                    calculator=Gaussian(
                        mem=SUBSYSTEM_MEM,
                        chk=CHECKPOINT_FILE,
                        save=None,
                        method=METHOD,
                        basis=BASIS,
                        scf="qc"
                    )))
            all_dft_traj+=dft_traj

            pm6_traj = self \
                .generate_subsets(atoms, Potential(
                    calculator=MOPAC(
                        method="PM6"
                    )))

            all_pm6_traj+=pm6_traj

            atom_types+=atoms.get_chemical_symbols()

        atom_types = list(set(atom_types))            
            
        # Figure out some sort of way to merge
        
        self.dft_model_file = self.train_model(atom_types, all_dft_traj)
        self.pm6_model_file = self.train_model(atom_types, all_pm6_traj)
    
    def subfrag(self, atoms: Atoms):
        dir_name=self.pdb_id
        project_dir=self.data_dir
        frag_dir=self.subfrag_dir
        pdb_file="/tmp/" + dir_name + ".pdb"
        com_file="/tmp/" + dir_name + ".com"
        gjf_file="/tmp/" + dir_name + ".gjf"
        sh_file="/tmp/" + dir_name + ".sh"
        write(pdb_file, atoms)
        os.write(sh_file, "#!/bin/bash\nmodule use /work2/01114/jfonner/frontera/modulefiles;  module load gaussian ;newzmat " + " -ipdb -ocom " + pdb_file + " " + com_file)
        os.system("chmod +x" + sh_file)
        os.system(sh_file)
        prepend = "%chk=" + dir_name + "\n" + "%nproc=56\n%njobs=6\n%Gver=g16\n%mem=5gb\n# pm6 \ngebf{dis=3, maxsubfrag=11, frag=protein}\n"
        content = os.read(com_file)
        append = "\n"
        os.write(gjf_file, prepend + content + append)
        
        os.system("lsqc " + gjf_file)
        os.write(sh_file, content)
        os.system("chmod +x " + sh_file)
        os.system(sh_file)
        
        with open(project_dir + "/" + dir_name + ".frg", "r") as file:
            index = 0
            for line in file:
                subsys = self.parse_fragment(line, atoms)
                write(frag_dir + "/" +  uuid.uuid1() + "_" + index + ".xyz", subsys)
                index+=1

        return frag_dir

    def parse_fragment(self, line: str, system: Atoms, frag_dir: str) -> Atoms:
        fields = line.split()
        index = fields[0]
        unknown_a = fields[1]
        subfrags_as_str = fields[2]
        unknown_b = fields[3]
        
        subfrags_as_array = subfrags_as_str \
            .replace(")", "") \
            .replace("(", "") \
            .split(",") \
        
        atoms_symbol = system.get_chemical_symbol()
        atoms_pos = system.get_positions()
        atoms = []

        for subfrag in subfrags_as_array:
            start = subfrag.split("-")[0] - 1
            end = subfrag.split("-")[1] - 1
            _range = start - end
            for index in _range:
                atom = Atom(symbol=atoms_symbol[i
    @abc.abstractmethodndex + start], position=atoms_pos[index + start])
                atoms.append(atom)
                
        mol = Atoms(atoms)
        return mol

    def generate_subsets(self, atoms: Atoms, potential: Potential):      
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

    def train_model(self, model_file: str, atypes: list, traj: list):
        raise NotImplementedError("train_model has not been implemented.")
        