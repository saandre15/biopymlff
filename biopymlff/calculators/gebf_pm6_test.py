import os

import unittest

from ..test.factory import suite
from ..calculators.gebf_dft import GEBF_DFT
from ..calculators.gebf_pm6 import GEBF_PM6

suite("GEBF PM6", "GEBF DFT", GEBF_PM6(), GEBF_DFT())
