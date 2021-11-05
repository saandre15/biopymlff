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

from biopymlff.calculators.GEBF_ML import GEBF_ML


# @ref https://libatoms.github.io/GAP/gap_fitting_tutorial.html#train-our-GAP_3b-model-from-the-command-line
class GEBF_GAP(GEBF_ML):

    _deprecated=object()

    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.', pdb_id=None,
                 **kwargs):

        super().__init__(
            restart=restart, 
            ignore_bad_restart_file=ignore_bad_restart_file, 
            label=label, 
            atoms=atoms, 
            directory=directory, 
            pdb_id=pdb_id,
            ext_type="gap")

    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):

        dataset_dir = self.data_dir + "/gap_dataset"
        os.mkdir(dataset_dir)

        write(dataset_dir + '/train.xyz', traj[0::2])
        write(dataset_dir + '/validate.xyz', traj[1::2])

        
        os.system("""
        gap_fit atoms_filename= # input data in extended XYZ format
            gap={                              # start of descriptor and kernel spec
                soap                              # first descriptor is a SOAP
                    l_max=6 n_max=12                  # number of angular and radial basis functions for SOAP
                    atom_sigma=0.5                    # Gaussian smearing width of atom density for SOAP, in Angstrom
                    cutoff=5.0                        # distance cutoff in the kernel, in Angstrom
                    radial_scaling=-0.5               # exponent of atom density scaling, power of distance
                    cutoff_transition_width=1.0       # distance across which kernel is smoothly taken to zero, in Angstrom
                    central_weight=1.0                # relative weight of central atom in atom density for SOAP
                    n_sparse=8000                     # number of representative points, M in Sec. II
                    delta=0.2                         # scaling of kernel, per descriptor, here for SOAP it is per atom, in eV
                    covariance_type=dot_product       # form of kernel
                    zeta=4                            # power kernel is raised to - together with dot_product gives a polynomial kernel
                    sparse_method=cur_points          # choice of representative points, here CUR decomposition of descriptor matrix
            }                                 # end of descriptor and kernel spec
            default_sigma={0.002 0.2 0.2 0.0}  # default regularisation corresponding to energy, force, virial, hessian
            config_type_sigma={                # start of per configuration-group regularisation spec, using groups defined in the input data file
                isolated_atom:0.0001:0.01:1.0:0.0:
                rss_rnd:0.03:0.4:0.5:0.0:
                rss_005:0.02:0.3:0.4:0.0:
                rss_200:0.01:0.2:0.2:0.0:
                rss_3c:0.005:0.1:0.1:0.00:
                cryst_dist:0.0003:0.03:0.05:0.00:
                cryst_dist_hp:0.005:0.1:0.1:0.0:
                liq_P4:0.003:0.3:0.5:0.0:
                liq_network:0.003:0.3:0.5:0.0:
                2D:0.001:0.03:0.05:0.0:
                ribbons:0.01:0.5:0.2:0.0
            }                                 # end of per configuration-group regularisation spec
            energy_parameter_name=energy       # name of the key in the input data file corresponding to the total energy
            force_parameter_name=forces        # name of the key in the input data file corresponding to the forces
            virial_parameter_name=virial       # name of the key in the input data file corresponding to the virial stress
            gp_file=gap.xml                    # name of output potential XML file
            sparse_jitter=1.0e-8               # extra diagonal regulariser
            do_copy_at_file=F                  # copy input data into potential XML file?
            sparse_separate_file=T             # write representative point data into a separate file not in the main potential XML
            core_param_file=P_r6_innercut.xml  # name of XML file containing the baseline potential (QUIP format)
            core_ip_args={IP Glue}             # initialisation string to call baseline potential
        """)     
        
        
        