import os
import uuid

from ase.atoms import Atoms
from ase.atom import Atom
from ase.geometry.analysis import Analysis
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.gaussian import Gaussian
from ase.calculators.amber import Amber

from shutil import which

from biopymlff.util.AtomGraph import AtomGraph

class GEBF(FileIOCalculator):
    
    implemented_properties = ['energy', 'forces']
    command = 'LSQC PREFIX.gjf'

    def __init__(self, restart=None,
                 ignore_bad_restart_file=Calculator._deprecated,
                 label=None, atoms=None,**kwargs):
        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, kawgs=kwargs)
        self.directory = os.getcwd() + "/data/" + label
        frg_file = self.get_fragment_file(atoms)
        os.rename(frg_file, self.directory + "/" + label + ".frg")

    def calculate(self, *args, **kwargs):
        lsqc = ('lsqc')
        if 'LSQC' in command:
            for program in lsqc:
                if which(program):
                    self.command = self.command.replace('LSQC', program)
                    break
            else: 
                raise EnvironmentError("lsqc is not installed on the system.")
        FileIOCalculator.calculate(self, args, kwargs)
    
    def subsystems(self, atoms: Atoms):
        exist = os.path.exists(self.directory + "/" + self.label + "_subsys")
        if not exist: raise IOError("Subsystem has not been created by LSQC yet.")
        

    def fragment(self, atoms: Atoms) -> list:
        G = AtomGraph()
        fragments = G.fragments_by_bond('C', 'C')
        return fragments
        
    def get_fragment_file(self, atoms: Atoms) -> str:
        filename = "/tmp/" + uuid.uuid1().hex + ".frg"
        fragment_file = open(filename, "a")
        fragments = self.fragment(atoms)
        for fragment in fragments:
            mol: Atoms = fragments
            G = AtomGraph(mol)
            spin_muliplicity = G.get_spin_multiplcity()
            charge = 0
            serial = 0
            fragment_indexes = []
            for atomI in mol:
                for atomJ in atoms:
                    if atomI == atomJ and atomI.symbol != 'H':
                        fragment_indexes.append(serial)
                    serial+=1
                serial = 0
            # Index Spin Muliplicity Fragment Indexes Charge
            line = f"""{serial} {spin_muliplicity} ({''.join(fragment_index).replace(' ', ',')}) 0"""
            fragment_file.write(line + "\n")
        fragment_file.close()
        return filename
                    
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
        
        with open(xyz_file, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(xyz_file, 'w') as fout:
            fout.writelines(data[2:])
        
        # Converts the pdb to a gaussian input file
        script = open(mk_gassuian_input, "w")
        script.write(f"""#!/bin/bash
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian
newzmat -ixyz -ocom {xyz_file} {com_file}"""
        )
        script.close()
        os.system("chmod +x " + mk_gassuian_input)
        os.system(mk_gassuian_input)
        # Converts a gaussian input file to a lsqc input file
        com_file_handler=open(com_file, "rt")
        com_file_content = com_file_handler.read()
        com_file_handler.close()
        script = open(gjf_file, "w")
        script.write(f"""%chk={dir_name}.chk
%nproc={cpu_count}
%njobs=6
%Gver=g16
%mem=10gb
# pm7

gebf{{dis=3 maxsubfrag=11 frag=protein}}

{com_file_content}

""")
        script.close()

        # Moves the input file and runs it to have the dataset within the project directory
        # xyz_file="/tmp/" + dir_name + ".xyz"
        script = open(run_lsqc, "w")
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
        script.close()

        os.system("chmod +x " + run_lsqc)
        os.system(run_lsqc)
        
        os.system("cd " + os.path.dirname(project_dir) + "; echo $PWD ;lsqc " + os.path.basename(gjf_file))

    def read_results(self):
        pass
            