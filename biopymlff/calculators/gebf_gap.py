import os
import sys
import uuid
import random

import numpy as np

import toml as toml
import os as os

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

from biopymlff.calculators.gebf_ml import GEBF_ML
from biopymlff.util.getenv import getenv
from biopymlff.util.gaussian import get_gaussian


class GEBF_GAP(GEBF_ML):
    """
    Generalized Based Energy Fragmentation utilizing Gaussian Approximate Potential 

    Notes
    -----
    dasdasasda

    References
    ----------
    https://libatoms.github.io/GAP/gap_fitting_tutorial.html#train-our-GAP_3b-model-from-the-command-line

    Bartók, Albert P., and Gábor Csányi. “Gaussian Approximation Potentials: A Brief Tutorial Introduction.” International Journal of Quantum Chemistry, vol. 115, no. 16, 2015, pp. 1051–1057., https://doi.org/10.1002/qua.24927. 
    """

    _deprecated=object()

    def __init__(self, restart=None, ignore_bad_restart_file=_deprecated,
                 label=None, atoms=None, directory='.', library=None
                 **kwargs):
        descriptor = SOAP(
            r_cutoff=gap_params["soap_r_c"],
            atom_sigma=gap_params["atom_sigma"],
            zeta=gap_params["zeta"],
            l_max=gap_params["l_max"],
            n_max=gap_params["N_R_l"],
            radial_scaling=-0.5,
            cutoff_trans_width=1.0,
            central_weight=1.0,
            n_sparse=8000,
            delta=0.2,
            covariance_type="dot_product",
            sparse_method="cur_points"
        )
        super(GEBF_ML, self).__init__(
            descriptors=[descriptor],
            restart=restart, 
            ignore_bad_restart_file=ignore_bad_restart_file, 
            label=label, 
            atoms=atoms, 
            directory=directory,
            library=library,
            ext_type="gap")


    def train_model(self, model_file: str, atypes: list, traj: list, type="default"):
        
        dataset_dir = self.directory + "/"

        xyz_file="/tmp/" + uuid.uuid1().hex + ".xyz"

        for a in traj:
            atoms: Atoms = a
            if type == "pm6":
                atoms.calc = get_gaussian(self.label, gaussian_params['pm6_method'])      
            elif type == "dft":
                atoms.calc = get_gaussian(self.label, gaussian_params['dft_method'], gaussian_params['dft_basis'])

        write(xyz_file, traj)

        soap_params = getenv()['gap']

        descriptors = ""

        # Make sure this functionatily works 
        counter = 0
        size = len(self.descriptors)
        for desc in self.descriptors:
            descriptors += str(desc) + (":" if counter != size else "")
        
        os.system(f"""
        gap_fit atoms_filename={xyz_file} # input data in extended XYZ format
            gap={{                              # start of descriptor and kernel spec
                {descriptors}
            }}                                 # end of descriptor and kernel spec
            default_sigma={{0.002 0.2 0.2 0.0}}  # default regularisation corresponding to energy, force, virial, hessian
            config_type_sigma={{                # start of per configuration-group regularisation spec, using groups defined in the input data file
                # isolated_atom:0.0001:0.01:1.0:0.0:
                # rss_rnd:0.03:0.4:0.5:0.0:
                # rss_005:0.02:0.3:0.4:0.0:
                # rss_200:0.01:0.2:0.2:0.0:
                # rss_3c:0.005:0.1:0.1:0.00:
                # cryst_dist:0.0003:0.03:0.05:0.00:
                # cryst_dist_hp:0.005:0.1:0.1:0.0:
                # liq_P4:0.003:0.3:0.5:0.0:
                # liq_network:0.003:0.3:0.5:0.0:
                # 2D:0.001:0.03:0.05:0.0:
                # ribbons:0.01:0.5:0.2:0.0
            }}                                 # end of per configuration-group regularisation spec
            energy_parameter_name=energy       # name of the key in the input data file corresponding to the total energy
            force_parameter_name=force        # name of the key in the input data file corresponding to the forces
            virial_parameter_name=virial       # name of the key in the input data file corresponding to the virial stress
            gp_file={model_file}                    # name of output potential XML file
            sparse_jitter=1.0e-8               # extra diagonal regulariser
            do_copy_at_file=F                  # copy input data into potential XML file?
            sparse_separate_file=T             # write representative point data into a separate file not in the main potential XML
            core_param_file=P_r6_innercut.xml  # name of XML file containing the baseline potential (QUIP format)
            core_ip_args={{IP Glue}}             # initialisation string to call baseline potential
        """)     
        
        
        