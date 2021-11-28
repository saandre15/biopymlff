import os

import unittest

from ..test.general import General_Test
from ..train.library import ProteinLibrary
from ..calculators.gebf_dft import GEBF_DFT
from ..calculators.gebf_pm6 import GEBF_PM6

class GEBF_PM6_Test(General_Test):

    
    def __init__(self, methodName: str):
        super().__init__(
            methodName,
            os.getcwd() + "/data/4znn.pdb",
            "GEBF PM6",
            "GEBF DFT",
            GEBF_PM6(),
            GEBF_DFT()
        )

if __name__ == '__main__':
    unittest.main()