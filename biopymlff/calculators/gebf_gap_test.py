import os

from biopymlff.test.calculator import Calculator_Test
from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_gap import GEBF_GAP
from biopymlff.test.factory import suite

class GEBF_GAP_Test(General_Test):
        
    def setUp(self):
        self.init_variables(os.path.join(proteins_path, "4znn.pdb"), "GEBF PM6", "GEBF DFT", GEBF_GAP(library=ProteinLibrary()), GEBF_DFT())
        super().setUp()