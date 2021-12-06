import os

from biopymlff.test.calculator import Calculator_Test
from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_gap import GEBF_GAP
from biopymlff.test.factory import suite

suite("GEBF GAP", "GEBF DFT", GEBF_GAP(library=ProteinLibrary()), GEBF_DFT())