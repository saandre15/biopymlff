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

from biopymlff.data.AtomGraph import AtomGraph


class GEBF(FileIOCalculator):
    """
    Generalized Energy Based Fragmentation
    """
    
    implemented_properties = ['energy', 'forces']
    command = 'LSQC PREFIX.gjf'

    def __init__(self, restart=None,
                 ignore_bad_restart_file=Calculator._deprecated,
                 label=None, atoms=None,**kwargs):
        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, kawgs=kwargs)
        self.directory = os.getcwd() + "/data/" + label
        if not os.path.exists(self.directory): os.mkdir(self.directory)
        self.frg_file = self.get_fragment_file(atoms)
        shutil.move(self.frg_file, self.directory + "/" + label + ".frg")

    def read_energy(self, lso_filepath='lso', gebf_filepath='[a]gebf', subsystems_dir="_subsys"):
        try:
            with open(lso_filepath) as file:
                lines = file.readlines()
                index_start = -1
                index_end = -1
                counter = 0
                for line in lines:
                    if "# GEBF-GETSUBS" in line:
                        index_start = counter
                    if "==========================================" in line and index_start != -1:
                        index_end = counter
                    counter+=1
                
                subsys_info=lines[index_start:index_end]
                index_start = -1
                index_mid = -1
                index_end = -1
                counter = 0
                for line in subsys_info:
                    if "Primitive subsystems":
                        index_start = counter
                    if "Derivative subsystems":
                        index_mid = counter
                    if "=======================================":
                        index_end = counter
                    counter+=1
                primative_subsystems = subsys_info[index_start + 1 : index_mid - 1]
                derivaitive_subsystem = subsys_info[index_mid + 1 : index_end - 1]
                subsystems = primative_subsystems + derivaitive_subsystem
                coefficents = []
                for line in subsystems:
                    vals = line.split()
                    coefficent = vals[1]
                    coefficents.append(coefficent)
                try:
                    with open(gebf_filepath) as file:
                        subsys_energy = []
                        lines: str = file.read()
                        subsystems: list = lines.split("# GEBF Subsys")
                        for subsys in subsystems:
                            subsys: str = subsys
                            subsys = subsys.split("# END")[0]
                            for line in subsys.split("\n"):
                                if self.which_potential() in line:
                                    energy= line.replace(self.which_potential(), '').strip()
                                    subsys_energy.append(energy)
                        
                        subsys_atoms = []
                        try:
                            subsys_files = os.listdir(subsystems_dir)
                            subsys_files = filter(lambda name: name.find(".sxyz") != -1, subsys_files)
                            for file in subsys_files:
                                subsys_atom = read_xyz(file, 0)
                                subsys_atoms.append(subsys_atom)
                            self.results['energy'] = self.calculate_potential_energy(coefficents, subsys_energy, subsys_atoms, atoms)
                        except IOError:
                            print("Failed to find subsystem directory containing topology information.")
                            raise ReadError()
                except IOError:
                    print("Unable to read gebf file for subsystem potential energy calculation.")
                    raise ReadError()
        except IOError:
            print("Unable to read lso file to find subsystem coefficents.")
            raise ReadError()
    
    def read_forces(self, force_filepath: str):
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

    def which_potential():
        return NotImplementedError("Define set potential to find. Example[E(SCF=)]")

    def calculate_potential_energy(self, coefficents: list, subsys_energy: list, subsys_atoms: list, atoms: Atoms):
        if len(coefficents) != len(subsys_energy) or len(subsys_energy) != len(atoms): raise IndexError("Make sure the length of each array and atom size are equal.")
        size=len(coefficents)
        total_subsys=0
        total_long_range=0
        for index in size:
            total_subsys+=coefficents[index]* (self.calculate_subsystem_pe(inital_subsys_energy, atoms))
        for a in atoms:
            for b in atoms:
                atomA: Atom = a
                atomB: Atom = b
                total_long_range += self.calculate_long_range_energy(a, b)
        return total_subsys + total_long_range
            
    def calculate_subsystem_pe(self, init_subsys_energy: float, atoms: Atoms):
        return init_subsys_energy
                       
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
        pass # call amber to get coefficent

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
        lso_filepath = os.getcwd() + "/" + self.label + "/" + self.label + "/" + self.label + ".lso"
        gebf_filepath = os.getcwd() + "/" +  self.label + "/" + self.label + "/" + self.label + ".gebf"
        subsystem_filepath = os.getcwd() + "/"  + self.label + "/" + self.label + "_subsys"
        force_filepath = os.getcwd() + "/" + self.label + "/" + self.label + "/" + self.label + ".force"
        self.read_energy(lso_filepath=lso_filepath, gebf_filepath=gebf_filepath, subsystems_dir=subsystem_filepath)
        self.read_forces(force_filepath=force_filepath)
            