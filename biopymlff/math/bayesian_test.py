import unittest
import os

from ase.atoms import Atoms 
from ase.io.proteindatabank import read_proteindatabank

from bayesian import get_alpha_beta, get_sigma

class Bayesian_Test(unittest.TestCase):
    
    def setUp(self):
        self.mol = read_proteindatabank(os.path.join(os.getcwd(), "data", "systems", "4znn.not_wat.pdb"))
        self.fingerprint = SOAP(l_max=6, n_max=3, atom_sigma=0.4, r_cutoff=3, radial_scaling=1, cutoff_trans_width=1, central_weight=1,
            n_sparse=1, delta=1, covariance_type="dot_product", zeta=2.5, species=mol.get_chemical_symbols(), sparse_method="cur_points")
        # Represents Fingerprints
        # NOTE: Algorithm begins to fail at around 50 atoms
        # TODO: Not working because matrix is not invertible
        self.x = self.fingerprint.to_tensor(self.mol)
        # Represents Potential Energy
        self.y = [random.randint(1, 30) for val in np.zeros(shape=(len(x), 1))]
    
    def test_alpha_beta(self):
        alpha, beta = get_alpha_beta(self.x, self.y)
        print("ALPHA " + alpha)
        print("BETA " + beta)

    def test_sigma(self):
        sigma = get_sigma(self.x, self.y)
        print("SIGMA " + sigma)

