import math

from ase.atoms import Atoms
from ase.calculators.calculator import Calculator

from ..util.getenv import getenv
from ..calculators.gebf import GEBF

class GEBF_PM6(GEBF):
    
    def __init__(self, restart=None, ignore_bad_restart_file=Calculator._deprecated, label=None, atoms=None, **kwargs):
        gaussian_params = getenv()["gaussian"]
        super(GEBF, self).__init__(
            method=gaussian_params["pm6_method"], 
            basis=None, 
            restart=restart, 
            ignore_bad_restart_file=ignore_bad_restart_file, 
            label=label, 
            atoms=atoms,
            kwargs=kwargs
        )
        

    