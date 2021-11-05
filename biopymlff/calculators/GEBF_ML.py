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


from biopymlff.calculators.ML import ML


class GEBF_ML(ML):

    implemented_properties = ['energy', 'energies', 'forces', 'stresses']

    _deprecated=object()

    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.', pdb_id=None, ext_type=None,
                 **kwargs):
        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, directory=directory, pdb_id=pdb_id)

        self.ext_type = ext_type

        self.dft_model_file=self.data_dir + "/dft_model.{self.ext_type}.xml"
        self.pm6_model_file=self.data_dir + "/pm6_model.{self.ext_type}.xml"
        self.add_model(self.dft_model_file)
        self.add_model(self.pm6_model_file)

    def get_subfrag_dir(self): return self.data_dir + "/" + self.pdb_id + "_subsys"

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
        # Runs some kind of structal optimization if neccesiary [MOPAC]
        self.subfrag(atoms)
        filename_woext=""
        numfiles = len([f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f)) and f[0] != '.'])
        
        frag_count = numfiles

        all_dft_traj = []
        all_pm6_traj = []
        atom_types = []

        for index in range(0, frag_count):
            atoms: Atoms = read_xyz(self.get_subfrag_dir() + filename_woext + "_" + index + ".xyz")
            
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
        print(dir(self))
        dir_name=self.pdb_id
        project_dir=self.data_dir
        frag_dir=self.get_subfrag_dir()
        xyz_file="/tmp/" + dir_name + ".xyz"
        com_file="/tmp/" + dir_name + ".com"
        gjf_file="/tmp/" + dir_name + ".gjf"
        
        mk_gassuian_input="/tmp/" + dir_name + ".gaussian.sh"
        mk_lsqc_input="/tmp/" + dir_name + ".lsqc.sh"
        run_lsqc="/tmp/" + dir_name + ".run.sh"
        cpu_count=os.cpu_count()
        # Write the atom to a pdb
        write(xyz_file, atoms)
        
        with open(xyz_file, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(xyz_file, 'w') as fout:
            fout.writelines(data[2:])
        
        # Converts the pdb to a gaussian input file
        script = open(mk_gassuian_input, "w")
        script.write(f"""#!/bin/bash
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian
newzmat -icom -ocom {xyz_file} {com_file}"""
        )
        script.close()
        os.system("chmod +x " + mk_gassuian_input)
        os.system(mk_gassuian_input)
        # Converts a gaussian input file to a lsqc input file
        com_file_handler=open(com_file, "rt")
        com_file_content = com_file_handler.read()
        com_file_handler.close()
        script = open(mk_lsqc_input, "w")
        script.write(f"""%chk={dir_name}.chk
%nproc={cpu_count}
%njobs=6
%Gver=g16
%mem=10gb
# pm7

# gebf{{dis=3, maxsubfrag=11, frag=protein}}

{com_file_content}

""")
        script.close()
        os.system("chmod +x " + mk_lsqc_input)
        os.system(mk_lsqc_input)

        # Moves the input file and runs it to have the dataset within the project directory
        # xyz_file="/tmp/" + dir_name + ".xyz"
        script = open(run_lsqc, "w")
        script.write(f"""#!/bin/bash
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian
export OMP_NUM_THREADS={cpu_count}
cd {project_dir}
mkdir {dir_name}
cp {xyz_file} {dir_name}
cp {gjf_file} {project_dir}
lsqc {os.path.basename(gjf_file)}""")
        script.close()

        os.system("chmod +x " + run_lsqc)
        os.system(run_lsqc)
        
        # Creates the subsystems
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
                atom = Atom(symbol=atoms_symbol[index + start], position=atoms_pos[index + start])
                atoms.append(atom)
                
        mol = Atoms(atoms)
        return mol

    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):
        raise NotImplementedError("train_model has not been implemented.")
        