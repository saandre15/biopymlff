import os

from biopymlff.test.calculator import Calculator_Test
from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_dp import GEBF_DP
from biopymlff.test.factory import suite

suite("GEBF DP", "GEBF DFT", GEBF_DP(library=ProteinLibrary()), GEBF_DFT())