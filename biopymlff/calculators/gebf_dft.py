from ase.calculators.calculator import Calculator

from biopymlff.calculators.gebf import GEBF
from biopymlff.util.getenv import getenv

class GEBF_DFT(GEBF):

    def __init__(self, restart=None, ignore_bad_restart_file=None, label=None, atoms=None, directory=".", **kwargs):
        gaussian_params = getenv()["gaussian"]
        super(GEBF, self).__init__(
            method=gaussian_params["dft_method"],
            basis=gaussian_params["dft_basis"], 
            restart=restart,
            ignore_bad_restart_file=ignore_bad_restart_file,
            label=label,
            atoms=atoms,
            directory=directory,
            kwargs=kwargs
        )
    
    