from ase.atoms import Atoms 

from ..


# NOTE: Algorithm begins to fail at around 50 atoms
mol: Atoms = read_proteindatabank(os.getcwd() + "/data/systems/4znn.not_wat.pdb")
# Potential Kernel for each atom
# TODO: Not working because matrix is not invertible
fingerprint = SOAP(l_max=6, n_max=3, atom_sigma=0.4, r_cutoff=3, radial_scaling=1, cutoff_trans_width=1, central_weight=1,
    n_sparse=1, delta=1, covariance_type="dot_product", zeta=2.5, species=mol.get_chemical_symbols(), sparse_method="cur_points")
x = fingerprint.to_tensor(mol)
y = [random.randint(1, 30) for val in np.zeros(shape=(len(x), 1))]

sigma = get_sigma(x, y)
print("SIGMA " + str(sigma))