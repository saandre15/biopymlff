from biopymlff.calculators.gaussian_gpu import Gaussian_GPU
from biopymlff.util.getenv import getenv

class PM6_GPU(Gaussian_GPU):
    
    def __init__(self, label=None):
        
        gaussian_params = getenv()["gaussian"]
        super(Gaussian_GPU, self).__init__(
            label=label,
            xc=gaussian_params["pm6_method"], 
            basis=None,
        )
        