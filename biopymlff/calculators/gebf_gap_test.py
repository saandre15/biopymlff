import os

from biopymlff.test.calculator import Calculator_Test
from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_gap import GEBF_GAP

class GEBF_GAP_Test(Calculator_Test):
    
    def __init__():
        super(Calculator_Test, self).__init__(
            os.getcwd() + "/data/4znn.pdb",
            "GEBF GAP",
            "GEBF DFT",
            GEBF_GAP(library=ProteinLibrary()),
            GEBF_DFT()
        )