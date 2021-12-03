import os
import uuid
import math
import shutil
import subprocess
import multiprocessing

from ase.atoms import Atoms
from ase.atom import Atom

from ase.geometry.analysis import Analysis

from ase.calculators.calculator import FileIOCalculator, CalculatorError, ReadError, CalculationFailed, Calculator, all_changes
from ase.calculators.gaussian import Gaussian
from ase.calculators.amber import Amber

from ase.io.xyz import read_xyz
from ase.io import write

from shutil import which

from ..data.atom_graph import AtomGraph
from ..util.getenv import getenv


class GEBF(FileIOCalculator):
    """
    Generalized Energy Based Fragmentation
    """
    
    implemented_properties = ['energy', 'forces']
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
                                if "Mulliken Charges" in line: # NOTE: Temporary Solution. Need to figure out how to get NPA charges instead of muliken
                                    read_charge = True
                                    continue
                                elif "Cartesian Gradiant" in line:
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

    def get_gaussian(self, xc=None, basis=None):
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
            xc=method,
            basis=basis,
            scf='(noincfock,novaracc,fermi,maxcycle=3000,ndamp=64,xqc)',
            int='acc2e=12'
        )
    
    def subsystems(self):
        mypath=self.directory
        files = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
        files = filter(lambda file: ".xyz" in file or ".sxyz" in file, files)
        subsystems = []
        for file in files:
            mol = read_xyz(file, 0)
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
        
    def get_fragment_file(self, atoms: Atoms) -> str:
        filename = "/tmp/" + uuid.uuid1().hex + ".frg"
        try:
            with open(filename, "w") as fragment_file:
                fragments_as_indexes = self.fragments_as_indexes(atoms)
                fragments_as_atoms = self.fragments_as_atoms(atoms)
                # TODO: Index not matching
                serial = 1
                for index in range(0, len(fragments_as_indexes)):
                    atoms = fragments_as_atoms[index]
                    indexes = fragments_as_indexes[index]
                    spin_multiplicity = AtomGraph(atoms).get_spin_multiplicity()
                    charge = AtomGraph(atoms).get_charges()
                    # charge = "+1"
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
        shutil.copyfile(self.get_fragment_file(atoms), self.label + ".frg")
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
        
        write(self.label + '.gjf', atoms, properties=None,
              format='gaussian-in', parallel=False, **self.parameters)
        with open(self.label + ".gjf", "r") as file:
            content = file.read()
            content = content.replace("Gaussian input prepared by ASE", "gebf{{frag=read charge=NPA}}")
            content = content.replace("\nkwargs", "")
            cpu_count = multiprocessing.cpu_count()
            content = "%nproc=" + str(cpu_count) + "\n" + content
            content = "%njobs=10\n" + content
            temp_content = ""
            lines = content.split("\n")
            for line in lines:
                if "TV" not in line: temp_content += line + "\n"
            content = temp_content

        with open(self.label + ".gjf", "w") as file:
            file.write(content + "\n")
            # NOTE: Figure out how to solve the lsqc issue on TACC
            print("lsqc should be called")
            
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
                    print(vals)                       
                    elec_count = vals[2]
                    charge_count = vals[3]
                    print("Electron Count " + str(elec_count))
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
            print(overwrite)
            file.write(overwrite)

        shutil.copy(self.label + "/" + self.label + ".frg", self.directory)
        shutil.copy(self.label + "/" + self.label + ".gjf", self.directory)
        shutil.rmtree(self.label)
                        

    def read_results(self):
        
        labc_filepath = os.path.join(os.getcwd(), self.label, self.label, self.label + ".labc")
        gebf_parent=os.path.join(os.getcwd(), self.label, self.label + "_subsys")
        gebf_filepaths = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(gebf_parent, f))]
        gebf_filepaths = filter(lambda file: "gebf" in file or ".agebf" in file, files)
        force_filepath = os.path.join(os.getcwd(), self.label, self.label, self.label + ".force")
        self.read_energy(labc_filepath=labc_filepath, gebf_filepaths=gebf_filepaths)
        self.read_forces(force_filepath=force_filepath)
            