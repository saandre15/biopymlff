import os

from biopymlff.util.gaussian import get_gaussian
from biopymlff.test.protein import Protein_Test
from biopymlff.test.general import General_Test
from biopymlff.train.library import ProteinLibrary
from biopymlff.calculators.gebf import GEBF
from biopymlff.calculators.gebf_dft import GEBF_DFT
from biopymlff.calculators.gebf_dp import GEBF_DP
from biopymlff.test.factory import suite

suite("GEBF DFT", "DFT", GEBF_DFT(), get_gaussian(label))