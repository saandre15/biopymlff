from ase.atoms import Atmos

from biopymlff.calculators.GEBF import GEBF

class GEBF_PM6(GEBF):
    
    def __init__():
        GEBF.__init__()

    def calculate_subsystem_pe(self, inital_subsys_energy: float, subsys_atoms: Atoms):
        base = super().calculate_subsystem_pe(init_subsys_energy, subsys_atoms)
        long_range_energy = 0
        for a in subsys_atoms:
            for b in subsys_atoms:
                if a == b: break
                long_range_energy += self.calculate_long_range_energy(a, b)
        return base + long_range_energy
    
    def which_potential(): return "E(SCF)="
        