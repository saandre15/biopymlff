import math

from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.mopac import MOPAC

class GEBF_PM7(Calculator):

    def subsys(self, atoms: Atoms, maxsubfrags: int):
        # distance thershold
        zeta=3
        # Max Number of Fragments
        n=4

        fragments = super().subfrag(atoms)
        subsystems = []

        for fragA in fragments:
            for fragB in fragments:
                if fragA == fragB: continue
                dist = self.get_closest_atom_distance(fragA, fragB)
                if dist == zeta:
                    sys = self.merge_atoms(fragA, fragB)
                    # cap dangling bonds with hydrogen
                    subsystems.append(sys)

        # delete smaller subsystems that are included in large subsystems
        index_to_remove=[]
        for i in len(subsystems):
            for j in len(subsystems):
                if i == j: continue
                if len(subsystems[i].get_positions()) != len(subsystems[j].get_positions()): index_to_remove.append(i)
                
        # create derivative subsystems
        

        return subsystems, coefficents

    def get_closest_atom_distance(self, a: Atoms, b: Atoms):
        coordsA = a.get_positions()
        coordsB = b.get_positions()
        dist = math.inf
        for coordA in coordsA:
            for coordB in coordsB:
                rad = self.get_radius(coordA, coordB)
                if rad < dist: dist = rad


        return rad

    def calculate(self, atoms: Atoms):
        coefficents, subsystems = self.subsys(atoms)
        subsys_energy = 0
        long_range_energy = self.get_long_range_colombic_energy(atoms)
        
        # Perform Amber NPA charge calculations

        size = len(atoms.get_chemical_symbols())
        charges = atoms.get_charges()
        coord = atoms.get_positions()

        calc = MOPAC('PM7')
        
        for index in range(0, len(coefficents)):
            subsys = subsystems[index]
            subsys.calc = MOPAC()
            subsys_energy += coefficents[index] * (subsys.get_potential_energy() - self.get_long_range_colombic_energy(subsys))
        
        return subsys_energy + lr_energy

    def get_long_range_colombic_energy(self, atoms: Atoms):
        size = len(atoms.get_chemical_symbols())
        long_range_energy=0

        for i in range(0, size):
            for j in range(0, size):
                if i == j: continue
                radius = self.get_radius(coord[i], coord[j])
                d_co=self.get_pairwise_dispersion_coefficent(a, b)
                lr_energy = ((charges[i] * charges[j]) / (radius)) + ((d_co ^ 12) / (radius ^ 12)) - ((d_co ^ 6) / (radius ^ 6))
                long_range_energy+=lr_energy

        return long_range_energy

    def get_pairwise_dispersion_coefficent(self, a: str, b: str):
        pass
            