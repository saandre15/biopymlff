from dataclasses import dataclass

from dscribe.descriptors.soap import SOAP as DSOAP

from ase.atoms import Atoms

from biopymlff.descriptors.descriptor import Descriptor


@dataclass
class SOAP(Descriptor):
    """
    
    """
    l_max: int
    n_max: int
    atom_sigma: float
    r_cutoff: float
    radial_scaling: float
    cutoff_trans_width: float
    central_weight: float
    n_sparse: float
    delta: float 
    covariance_type: str
    zeta: int 
    sparse_method: str    
    species: list
    
    def __str__(self):
        return f"""soap
            l_max={self.l_max} n_max={self.n_max}
            atom_sigma={self.atom_sigma}
            cutoff={self.r_cutoff}
            radial_scaling={self.radial_scaling}
            cutoff_transitition_width={self.cutoff_trans_width}
            central_weight={self.central_weight}
            n_sparse={self.n_sparse}
            delta={self.delta}
            covariance_type={self.covariance_type}
            zeta={self.zeta}
            sparse_method={self.sparse_method}
        """

    def to_dscribe(self):
        return DSOAP(rcut=self.r_cutoff, nmax=self.n_max, lmax=self.l_max, sigma=self.atom_sigma, rbf="gto", weighting=None, 
            crossover=True, average='off', species=self.species, periodic=False, sparse=False, dtype='float64')

    def to_tensor(self, atoms: Atoms):
        return self.to_dscribe().create(atoms)

    def to_derivative_tensor(self, atoms: Atoms):
        return self.to_dscribe().derivatives(atoms)