import os

import unittest

from ..test.general import General_Test
from ..train.library import ProteinLibrary
from ..calculators.gebf_dft import GEBF_DFT
from ..calculators.gebf_pm6 import GEBF_PM6

class GEBF_PM6_Test(General_Test):

    @classmethod
    def setUpClass(cls):
        cls.init_variables(os.getcwd() + "/data/systems/4znn.pdb", "GEBF PM6", "GEBF DFT", GEBF_PM6(), GEBF_DFT())