import os

from ase.calculators.calculator import Calculator

from ..test.protein import Protein_Test
from ..test.lipid import Lipid_Test

def suite(source_method: str, target_method: str, source: Calculator, target: Calculator):

    path = os.path.join(os.getcwd(), "data", "systems")
    proteins_path = os.path.join(path, "proteins")
    lipids_path = os.path.join(path, "lipids")
    
    class Protein_Test_01(Protein_Test):
        
        def setUp(self):
            self.init_variables(os.path.join(proteins_path, "4znn.pdb"), source_method, target_method, source, target)
            super().setUp()

    class Protein_Test_02(Protein_Test):
        
        def setUp(self):
            self.init_variables(os.path.join(proteins_path, "1xq8.pdb"), source_method, target_method, source, target)
            super().setUp()

    class Protein_Test_03(Protein_Test):
        
        def setUp(self):
            self.init_variables(os.path.join(proteins_path, "2yvb.pdb"), source_method, target_method, source, target)
            super().setUp()
    
    class Lipid_Test_01(Lipid_Test):
        
        def setUp(self):
            self.init_variables(os.path.join(lipids_path, "4c7r.pdb"), source_method, target_method, source, target)
            super().setUp()

    
    return [Protein_Test_01, Protein_Test_02, Protein_Test_03, Lipid_Test_01]