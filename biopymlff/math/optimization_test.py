import unittest

import fingerprint

class Optimization_Test(unittest.TestCase):

    def setUp(self):
        self.mol = read_proteindatabank(os.getcwd() + "/data/systems/4znn.not_wat.pdb")
        # mol = Atoms(symbols='H2O')
        self.fingerprint = SOAP(l_max=6, n_max=3, atom_sigma=0.4, r_cutoff=3, radial_scaling=1, cutoff_trans_width=1, central_weight=1,
            n_sparse=1, delta=1, covariance_type="dot_product", zeta=2.5, species=mol.get_chemical_symbols(), sparse_method="cur_points").to_tensor(mol)
        self.x = self.fingerprint
        print(cur_approximation(fingerprint, np.linalg.matrix_rank(fingerprint)))

    def test_cur_matrix_approximation(self):
        