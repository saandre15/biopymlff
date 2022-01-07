from biopymlff.test.general import General_Test
from biopymlff.calculators.pm6_gpu import PM6_GPU
from biopymlff.calculators.dft_gpu import DFT_GPU

class PM6_GPU_Test(General_Test):
        
    def setUp(self):
        self.init_variables(os.path.join(proteins_path, "4znn.pdb"), "PM6 GPU", "DFT GPU", PM6_GPU(), DFT_GPU())
        super().setUp()