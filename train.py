import os

import ase

from biopymlff.calculators.GEBF_GAP import GEBF_GAP
# Complete 

pdb_id="4znn"
pdb_file = os.getcwd() + "/data/systems/" + pdb_id + ".pdb"

atoms = ase.io.read(pdb_file)

calc = GEBF_GAP(atoms=atoms, pdb_id=pdb_id)

calc.train()





