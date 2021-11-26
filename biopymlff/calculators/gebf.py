import os
import uuid
import math
import shutil

from ase.atoms import Atoms
from ase.atom import Atom

from ase.geometry.analysis import Analysis

from ase.calculators.calculator import FileIOCalculator, CalculatorError, ReadError, CalculationFailed, Calculator
from ase.calculators.gaussian import Gaussian
from ase.calculators.amber import Amber

from ase.io.xyz import read_xyz

from shutil import which

from biopymlff.data.atom_graph import AtomGraph


class GEBF(FileIOCalculator):
    """
    Generalized Energy Based Fragmentation
    """
    
    implemented_properties = ['energy', 'forces']
    command = 'LSQC PREFIX.gjf'

    def __init__(self, method: str, basis: str=None, restart=None,
                 ignore_bad_restart_file=Calculator._deprecated,
                 label=None, atoms=None, directory='.',**kwargs):
        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, directory=directory, kawgs=kwargs)
        self.directory = os.getcwd() + "/data/" + label
        if not os.path.exists(self.directory): os.mkdir(self.directory)
        self.frg_file = self.get_fragment_file(atoms)
        shutil.move(self.frg_file, self.directory + "/" + label + ".frg")

    def read_energy(self, labc_filepath='labc', gebf_filepaths: list=[]):
        if len(gebf_filepaths) != len(xyz_filepaths):
            raise ReadError("Cannot obtain potential energy where gebf file size doesn't match xyz file size.")
        try:
            with open(labc_filepath) as file:
                lines = file.readlines()
                read_coefficient = False
                coefficents = []
                for line in lines:
                    if "coef:" in line:
                        read_coefficient = True
                        continue
                    if "Table: Lab for subsystems" in line:
                        read_coefficient = False 
                        break
                    coefs = line.split()
                    for coef in coefs:
                        coefficents.append(coef)

                try:
                    subsys_atoms = self.subsystems()    
                    counter = 0
                    for path in gebf_filepaths:
                        with open(gebf_filepaths) as file:
                            lines = file.readlines()
                            atoms: Atoms = subsys_atoms[counter]
                            read_charge = False
                            charges = []
                            for line in lines:
                                if "NPA Charges" in line: 
                                    read_charge = True
                                    continue
                                if "Cartesian Gradiant" in line:
                                    read_charge = False
                                    break
                                crgs = line.split()
                                for charge in crgs:
                                    charges.append(charge)
                            counter+=1
                            atoms.set_initial_charges(charges)
                    self.calculate_potential_energy(coefficents, subsys_atoms, atoms)
                except IOError:
                    print("Unable to read gebf file for subsystem potential energy calculation.")
                    raise ReadError()
        except IOError:
            print("Unable to read lso file to find subsystem coefficents.")
            raise ReadError()
    
    def read_forces(self, force_filepath: str=None):
        try:
            with open(force_filepath) as file:
                self.results['forces'] = []
                forces = file.readlines()
                f = forces.split()
                f_x = f[0]
                f_y = f[1]
                f_z = f[2]
                f = (f_x, f_y, f_z)
                self.results['forces'].append(f)
        except IOError:
            print("Unable to read force file. Make sure the force_filepath is correct.")
            raise ReadError()

    def calculate_potential_energy(self, coefficents: list, subsys_atoms: list, atoms: Atoms):
        if len(coefficents) != len(subsys_atoms): raise IndexError("Make sure the length of each array and atom size are equal.")
        size=len(coefficents)
        total_subsys=0
        total_long_range=0
        for index in size:
            total_subsys+=coefficents[index]* (self.calculate_subsystem_pe(0, subsys_atoms[index]))
        for a in atoms:
            for b in atoms:
                atomA: Atom = a
                atomB: Atom = b
                total_long_range += self.calculate_long_range_energy(a, b)
        return total_subsys + total_long_range
            
    def calculate_subsystem_pe(self, init_subsys_energy: float, subsys_atoms: Atoms):
        gap_params = getenv()["gap"]
        cutoff_radius = gap_params["soap_r_c"]
        subsys_energy = init_subsys_energy
        for a in atoms:
            for b in atoms:
                if a == b: continue
                f_cutoff = 1 if cutoff_radius < radius else 0.5 * ( 1 - math.cos(math.pi * radius * (1 / cutoff_radius)) )
                val = self.calculate_long_range_energy(a, b) * f_cutoff
                subsys_energy+=val
        return subsys_energy
                       
    def calculate_long_range_energy(self, a: Atom, b: Atom):
        charge_energy = 0
        van_der_val_energy = 0
        charge_energy += (a.charge * b.charge) / self.get_radius(a, b)
        van_der_val_energy += (self.get_lorentz_berthetot_coefficent(a, b) ^ 12 / self.get_radius(a, b) ^ 12) \
            - (self.get_lorentz_berthetot_coefficent(a, b) ^ 12 / self.get_radius(a, b) ^ 6)
        return charge_energy + van_der_val_energy

    def get_radius(self, a: Atom, b: Atom):
        return math.sqrt((a.x - b.x) ^ 2 + (a.y - b.y) ^ 2 + (a.z + b.z) ^ 2)

    def get_lorentz_berthetot_coefficent(self, a: Atom, b: Atom):
        # TODO figure out how to read the param file and perform coeffiecne tcalculation
        amber_params = getenv()["amber"]
        return

    def calculate(self, atoms: Atoms):
        lsqc = ('lsqc')
        if 'LSQC' in self.command:
            for program in lsqc:
                if which(program):
                    print("test")
                    self.command = self.command.replace('LSQC', program)
                    break
            else: 
                raise EnvironmentError("lsqc is not installed on the system.")
        FileIOCalculator.calculate(self, args, kwargs)
    
    def subsystems(self):
        subsys_dir = self.directory + "/" + self.label + "_subsys"
        exist = os.path.exists(subsys_dir)
        if not exist: raise IOError("Subsystem has not been created by LSQC yet.")
        files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
        size = len(files) / 4
        subsystems = []
        for index in len(1, size + 1):
            name = self.label + "_" + index
            # gaus = Gaussian(label=name, directory=self.directory + "/" + self.label + "_subsys")
            mol: Atoms = read_xyz(name + ".sxyz", 0)
            subsystems.append(mol)
        return subsystems

    def fragments_as_indexes(self, atoms: Atoms) -> list:
        G = AtomGraph(atoms)
        fragments = G.fragments_by_bond_as_indexes('C', 'C')
        return fragments
    
    def fragments_as_atoms(self, atoms: Atoms) -> list:
        G = AtomGraph(atoms)
        fragments = G.fragment_by_bond_as_atoms_list('C', 'C')
        return fragments
    
    def get_spin_multiplcity(self, atoms: Atoms):
        return AtomGraph(atoms).get_spin_multiplicity()
        
    def get_fragment_file(self, atoms: Atoms) -> str:
        filename = "/tmp/" + uuid.uuid1().hex + ".frg"
        # print(atoms)
        try:
            with open(filename, "a") as fragment_file:
                # print(atoms)
                fragments_as_indexes = self.fragments_as_indexes(atoms)
                print(fragments_as_indexes)
                fragments_as_atoms = self.fragments_as_atoms(atoms)
                serial = 1
                for index in range(0, len(fragments_as_indexes)):
                    atoms = fragments_as_atoms[index]
                    indexes = fragments_as_indexes[index]
                    spin_multiplicity = self.get_spin_multiplcity(atoms)
                    charge = 0
                    temp_indexes = []
                    for index in indexes:
                        temp_indexes.append(index + 1)
                    # Index Spin Muliplicity Fragment Indexes Charge
                    line = f"""{str(serial)} {str(spin_multiplicity)} ({','.join([str(val) for val in temp_indexes])}) 0"""
                    fragment_file.write(line + "\n")
                    serial+=1
                fragment_file.close()
                return filename
        except IOError: print("Unable to create a fragment file in GEBF.")
                    
    def write_input(self, atoms: Atoms):
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
        
        try:
            with open(xyz_file, 'r') as fin:
                data = fin.read().splitlines(True)
            with open(xyz_file, 'w') as fout:
                fout.writelines(data[2:])
        except IOError: 
            print("Failed to read XYZ file to write input.")
            raise ReadError()

        # Converts the pdb to a gaussian input file
        with open(mk_gassuian_input, "w") as script:
            script.write(f"""#!/bin/bash
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian
newzmat -ixyz -ocom {xyz_file} {com_file}"""
            )
            script.close()
            os.system("chmod +x " + mk_gassuian_input)
            os.system(mk_gassuian_input)
            # Converts a gaussian input file to a lsqc input file
            with open(com_file, "rt") as com_file_handler:
                com_file_content = com_file_handler.read()
                with open(gjf_file, "") as script:
                    
                    script.write(f"""%chk={dir_name}.chk
%nproc={cpu_count}
%njobs=6
%Gver=g16
%mem=10gb
# pm6

gebf{{frag=read}}

{com_file_content}

""")

                    # Moves the input file and runs it to have the dataset within the project directory
                    # xyz_file="/tmp/" + dir_name + ".xyz"
                    with open(run_lsqc, "w") as script:
                        script.write(f"""#!/bin/bash
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian
export OMP_NUM_THREADS={cpu_count}
cd {os.path.dirname(project_dir)}
mkdir {dir_name}
cp {xyz_file} {dir_name}
cp {gjf_file} .
mkdir -p {self.get_subfrag_dir()}
""")    

                        os.system("chmod +x " + run_lsqc)
                        os.system(run_lsqc)
                        
                        os.system("cd " + os.path.dirname(project_dir) + "; echo $PWD ;lsqc " + os.path.basename(gjf_file))

    def read_results(self):
        lso_filepath = os.getcwd() + "/" + self.label + "/" + self.label + ".labc"
        gebf_filepath = os.getcwd() + "/" +  self.label + "/" + self.label + ".gebf"
        force_filepath = os.getcwd() + "/" + self.label + "/" + self.label + ".force"
        self.read_energy(lso_filepath=lso_filepath, gebf_filepath=gebf_filepath)
        self.read_forces(force_filepath=force_filepath)
            