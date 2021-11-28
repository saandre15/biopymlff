import os

from biopymlff.test.calculator import Calculator_Test
from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_pm6 import GEBF_PM6

class GEBF_PM6_Test(Calculator_Test):
    
    def __init__():
        super(Calculator_Test, self).__init__(
            os.getcwd() + "/data/4znn.pdb",
            "GEBF PM6",
            "GEBF DFT",
            GEBF_PM6(),
            GEBF_DFT()
        )

if __name__ == '__main__':
    unittest.main()