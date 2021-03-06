import os
import uuid
import math
import shutil
import subprocess
import multiprocessing

import numpy as np

from ase.atoms import Atoms
from ase.atom import Atom

from ase.geometry.analysis import Analysis

from ase.calculators.calculator import FileIOCalculator, CalculatorError, ReadError, CalculationFailed, Calculator, all_changes
from ase.calculators.gaussian import Gaussian
from ase.calculators.amber import Amber

from ase.md.langevin import Langevin

from ase.io.xyz import read_xyz
from ase.io.proteindatabank import read_proteindatabank
from ase.io import write

from shutil import which

from ..data.atom_graph import AtomGraph
from ..util.getenv import getenv
from ..data.atom_graph_edge_type import AtomGraphEdgeType


class GEBF(FileIOCalculator):
    """
    Generalized Energy Based Fragmentation
    """
    
    implemented_properties = ['energy', 'forces', 'free_energy', 'energies']
    command = 'LSQC PREFIX.gjf'

    def __init__(self, restart=None,
                 ignore_bad_restart_file=Calculator._deprecated,
                 label=None, atoms=None, directory='.',**kwargs):
        """ 
        Takes ASE Gaussian Parameters As Well
        """
        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, directory=directory, kawgs=kwargs)
        if not os.path.exists(self.directory): os.mkdir(self.directory)
        self.frg_file = self.get_fragment_file(atoms)
        shutil.move(self.frg_file, self.directory + "/" + label + ".frg")


    def read_energy(self, labc_filepath='.labc', gebf_filepaths: list=[], charge_filepath: str='.cha'):
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
                    gebf_filepaths=sorted(gebf_filepaths)
                    for path in gebf_filepaths: # Sets subsystem charges
                        with open(path) as file:
                            lines = file.readlines()
                            atoms: Atoms = subsys_atoms[counter]
                            read_charge = False
                            charges = []
                            for line in lines:
                                #  if "NPA Charges" in line:
                                if "Mulliken Charges" in line: # NOTE: Temporary Solution. Need to figure out how to get NPA charges instead of muliken
                                    read_charge = True
                                    continue
                                elif "Cartesian Gradient" in line:
                                    read_charge = False
                                    break
                                if read_charge == True: 
                                    crgs = line.split()
                                    for charge in crgs:
                                        charges.append(float(charge))
                            counter+=1
                            atoms.set_initial_charges(charges)
                    with open(charge_filepath) as file:
                        lines = file.readlines()
                        charges = []
                        for line in lines:
                            vals = line.split(" ")
                            print(vals)
                            col = []
                            for val in vals:
                                if val != '':
                                    col.append(val)
                            print(val)
                            charge_val = float(col[2])
                            charges.append(charge_val)
                        self.atoms.set_initial_charges(charges)
                    energy = self.calculate_potential_energy(coefficents, subsys_atoms, self.atoms)
                    self.results["energy"]=float(energy)
                    self.results["free_energy"]=float(energy)
                    self.results["energies"] = np.empty(shape=(len(self.atoms), 1))
                    index = 0
                    for atom in self.atoms:
                        self.results["energies"][index] = self.calculate_interatomic_pe(coefficents, subsys_atoms, self.atoms, atom)
                        index+=1
                    # NOTE: If method above does not work. Then adding stress
                except IOError:
                    print("Unable to read gebf file for subsystem potential energy calculation or charge file for system charge.")
                    raise ReadError()
        except IOError:
            print("Unable to read lso file to find subsystem coefficents.")
            raise ReadError()
    
    def read_forces(self, force_filepath: str):
        try:
            with open(force_filepath) as file:
                self.results['forces'] = []
                lines = file.readlines()
                count = len(lines)
                index = 0
                self.results['forces'] = np.empty(shape=(count, 3))
                for line in lines:
                    f = line.split()
                    f_x = float(f[0].replace("D", "e"))
                    f_y = float(f[1].replace("D", "e"))
                    f_z = float(f[2].replace("D", "e"))
                    f = np.array([f_x, f_y, f_z])
                    self.results['forces'][index] = f
                    index+=1
                
        except IOError:
            print("Unable to read force file. Make sure the force_filepath is correct.")
            raise ReadError()

    def calculate_potential_energy(self, coefficents: list, subsys_atoms: list, atoms: Atoms):
        if len(coefficents) != len(subsys_atoms): raise IndexError("Make sure the length of each array and atom size are equal.")
        size=len(coefficents)
        total_subsys=0
        total_long_range=0
        for index in range(0, size):
            total_subsys+=float(coefficents[index])* float(self.calculate_subsystem_pe(0, subsys_atoms[index]))
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
        for a in subsys_atoms:
            for b in subsys_atoms:
                if a == b: continue
                radius = self.get_radius(a, b)
                f_cutoff = 1 if cutoff_radius < radius else 0.5 * ( 1 - math.cos(math.pi * radius * (1 / cutoff_radius)) )
                long_range_energy = self.calculate_long_range_energy(a, b)
                val = (long_range_energy * f_cutoff) - long_range_energy
                subsys_energy+=val
        return subsys_energy

    def calculate_interatomic_pe(self, coefficients: list, subsys_atoms: list, atoms: Atoms, target: Atom):
        if len(coefficients) != len(subsys_atoms): raise IndexError("Cannot calculate the interatomic potential energy where the coefficents length size does not match the subsystem length size.")
        gap_params = getenv()["gap"]
        total_subsys = 0
        total_long_range = 0
        cutoff_radius = gap_params["soap_r_c"]
        for index in range(0, len(coefficients)):
            coef = coefficients[index]
            subsys = subsys_atoms[index]
            subsys_energy = 0
            for atom in subsys:
                if atom == target: continue
                radius = self.get_radius(atom, target)
                f_cutoff = 1 if cutoff_radius < radius else 0.5 * ( 1 - math.cos(math.pi * radius * (1 / cutoff_radius)) )
                long_range_energy = self.calculate_long_range_energy(atom, target)
                val = (long_range_energy * f_cutoff) - long_range_energy
                subsys_energy += val
            total_subsys+=float(coefficients[index]) * float(subsys_energy)
        
        for atom in atoms:
            if atom == target: continue 
            total_long_range += self.calculate_long_range_energy(atom, target)

        return total_subsys + total_long_range  
                       
    def calculate_long_range_energy(self, a: Atom, b: Atom):
        charge_energy = 0
        van_der_val_energy = 0
        charge_energy += (a.charge * b.charge) / self.get_radius(a, b)
        van_der_val_energy += 0
        # van_der_val_energy += (self.get_lorentz_berthetot_coefficent(a, b) ** 12 / self.get_radius(a, b) ** 12) \
        #     - (self.get_lorentz_berthetot_coefficent(a, b) ** 6 / self.get_radius(a, b) ** 6)
        # return charge_energy + van_der_val_energy
        return charge_energy + van_der_val_energy

    def get_radius(self, a: Atom, b: Atom):
        return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z + b.z) ** 2)

    def get_lorentz_berthetot_coefficent(self, a: Atom, b: Atom):
        # TODO figure out how to read the param file and perform coeffiecne tcalculation
        # amber_params = getenv()["amber"]
        # amber_home = amber_params["amber_home"]
        return 1

    def calculate(self, *args, **kwargs):
        lsqc = ['lsqc']
        if 'LSQC' in self.command:
            for program in lsqc:
                if which(program):
                    self.command = self.command.replace('LSQC', program)
                    break
            else: 
                raise EnvironmentError("lsqc is not installed on the system.")
        FileIOCalculator.calculate(self, *args, **kwargs)

    def calculate_repair(self, atoms: Atoms, properties=[], system_changes=all_changes): 
        Calculator.calculate(self, atoms, properties, system_changes)
        # self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise CalculatorSetupError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        command = self.command
        if 'PREFIX' in command:
            command = command.replace('PREFIX', self.prefix)
        try:
            # NOTE: Not executing for what ever reason???
            proc = subprocess.Popen(command, shell=True, cwd=self.directory, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
        except OSError as err:
            # Actually this may never happen with shell=True, since
            # probably the shell launches successfully.  But we soon want
            # to allow calling the subprocess directly, and then this
            # distinction (failed to launch vs failed to run) is useful.
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err
        # stderr, stdout = proc.communicate()
        errorcode = proc.wait()

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                   '{} with error code {}'.format(self.name, command,
                                                  path, errorcode))
            raise CalculationFailed(msg)
    
    def subsystems(self):
        mypath=os.path.join(self.directory, self.label, self.label + "_subsys")
        files = [os.getcwd() + "/" + self.label + "/" + self.label + "_subsys/" + f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
        files = filter(lambda file: ".xyz" in file or ".sxyz" in file, files)
        files = sorted(files)
        subsystems = []
        for path in files:
            with open(path) as file:
                print(path)
                mol = Atoms(list(read_xyz(file, 0)))
                subsystems.append(mol)
        return subsystems

    def fragments_as_indexes(self, atoms: Atoms) -> list:
        print("testing")
        G = AtomGraph(atoms)
        fragments = G.fragments_by_bond_as_indexes('C', 'C', [AtomGraphEdgeType.SINGLE])
        return fragments
    
    def fragments_as_atoms(self, atoms: Atoms) -> list:
        G = AtomGraph(atoms)
        fragments = G.fragment_by_bond_as_atoms_list('C', 'C', [AtomGraphEdgeType.SINGLE])
        return fragments
        
    def get_fragment_file(self, atoms: Atoms) -> str:
        filename = os.path.join("/", "tmp", uuid.uuid1().hex + ".frg")
        if os.path.exists(filename): os.remove(filename)
        try:
            with open(filename, "w") as fragment_file:
                print("apple" + str(1))
                fragments_as_indexes = self.fragments_as_indexes(atoms)
                fragments_as_atoms = self.fragments_as_atoms(atoms)
                # TODO: Index not matching
                serial = 1
                for index in range(0, len(fragments_as_indexes)):
                    atoms = fragments_as_atoms[index]
                    indexes = fragments_as_indexes[index]
                    spin_multiplicity = AtomGraph(atoms).get_spin_multiplicity()
                    charge = AtomGraph(atoms).get_charges()
                    temp_indexes = []
                    for index in indexes:
                        temp_indexes.append(index + 1)
                    # Index Spin Muliplicity Fragment Indexes Charge
                    line = f"""{str(serial)} {str(spin_multiplicity)} ({','.join([str(val) for val in temp_indexes])}) {str(charge)}"""
                    fragment_file.write(line + "\n")
                    serial+=1
                fragment_file.close()
                return filename
        except IOError: print("Unable to create a fragment file in GEBF.")
                    
    def write_input(self, atoms: Atoms, properties=None, system_changes=None):
        print(atoms)
        print(self.get_fragment_file(atoms))
        shutil.copyfile(self.get_fragment_file(atoms), self.label + ".frg") # TODO: Fix missing self.label
        general_params = getenv()['general']
        gaussian_params = getenv()['gaussian']
        
        # TODO: Per subsystem memory not total memory
        # self.parameters["mem"] = gaussian_params['memory'] \
        #         if gaussian_params['memory'] != "auto" \
        #         else str(math.floor((os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024.**3)) - 5)) + "GB"
        self.parameters["mem"] = "3GB"
        self.parameters["chk"] = gaussian_params['checkpoint_file'] if gaussian_params['checkpoint_file'] != "auto" else self.label + ".chk"
        self.parameters["scf"]='noincfock,novaracc,fermi,maxcycle=3000,ndamp=64,xqc'
        self.parameters["int"]='acc2e=12'
        self.parameters["pop"] = gaussian_params['pop'] # NOTE: Testing natural population analysis for PM6 and DFT
        
        write(self.label + '.gjf', atoms, properties=None,
              format='gaussian-in', parallel=False, **self.parameters)
        with open(self.label + ".gjf", "r") as file:
            content = file.read()
            title = "gebf " + "{" + " frag=read charge=" + gaussian_params['pop'] + " }"
            content = content.replace("Gaussian input prepared by ASE", title)
            content = content.replace("\nkwargs", "")
            cpu_count = multiprocessing.cpu_count()
            content = "%nproc=" + str(cpu_count) + "\n" + content
            content = "%njobs=10\n" + content
            temp_content = ""
            lines = content.split("\n")
            for line in lines:
                if "TV" not in line: temp_content += line + "\n" # Removes periodic boundary condition
            content = temp_content

        with open(self.label + ".gjf", "w") as file:
            file.write(content + "\n")
            
        self.calculate_repair(atoms=atoms, properties=properties)
        with open(self.label + "/" + self.label + ".lso") as file:
            lines = file.readlines()
            prereading_mode = False
            reading_mode = False
            charge_correction = []
            for line in lines:
                if "============" in line and reading_mode == True:
                    print("this is called")
                    prereading_mode = False
                    reading_mode = False
                    break
                elif reading_mode == True:
                    vals = line.split()                  
                    elec_count = vals[2]
                    charge_count = vals[3]
                    if int(elec_count) % 2 == 1: 
                        charge_correction.append(int(charge_count)+1)
                    else: 
                        charge_correction.append(int(charge_count))
                elif "Frag NAtoms Elec Char Mult" in line:
                    prereading_mode = True
                elif "----------" in line and prereading_mode == True:
                    reading_mode = True 

        with open(self.label + "/" + self.label + ".frg", "r") as file:
            lines = file.readlines()
            overwrite = ""
            counter=0
            for line in lines:
                vals = line.split()
                overwrite+=vals[0] + " " + vals[1] + " " + vals[2] + " " + (str(charge_correction[counter]) if charge_correction[counter] < 1 else "+" + str(charge_correction[counter])) + "\n"
                counter+=1

        os.remove(self.label + "/" + self.label + ".frg")
        
        with open(self.label + "/" + self.label + ".frg", "w") as file:
            file.write(overwrite)

        shutil.copy(self.label + "/" + self.label + ".frg", self.directory)
        shutil.copy(self.label + "/" + self.label + ".gjf", self.directory)
        shutil.rmtree(self.label)
                        

    def read_results(self):
        labc_filepath = os.path.join(os.getcwd(), self.label, self.label, self.label + ".labc")
        gebf_parent=os.path.join(os.getcwd(), self.label, self.label + "_subsys")
        gebf_filepaths = [os.getcwd() + "/" + self.label + "/" + self.label + "_subsys/" + f for f in os.listdir(gebf_parent) if os.path.isfile(os.path.join(gebf_parent, f))]
        print(gebf_filepaths)
        gebf_filepaths = filter(lambda file: ".gebf" in file or ".agebf" in file, gebf_filepaths)
        force_filepath = os.path.join(os.getcwd(), self.label, self.label, self.label + ".force")
        charge_filepath = os.path.join(os.getcwd(), self.label, self.label, self.label + ".cha")
        self.read_energy(labc_filepath=labc_filepath, gebf_filepaths=gebf_filepaths, charge_filepath=charge_filepath)
        self.read_forces(force_filepath=force_filepath)
            