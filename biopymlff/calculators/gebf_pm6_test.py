import os

import unittest

from ..test.factory import suite
from ..test.protein import Protein_Test
from ..test.general import General_Test
from ..calculators.gebf_dft import GEBF_DFT
from ..calculators.gebf_pm6 import GEBF_PM6

path = os.path.join(os.getcwd(), "data", "systems")
proteins_path = os.path.join(path, "proteins")

class GEBF_PM6_Test(General_Test):
        
    def setUp(self):
        self.init_variables(os.path.join(proteins_path, "4znn.pdb"), "GEBF PM6", "GEBF DFT", GEBF_PM6(), GEBF_DFT())
        super().setUp()
