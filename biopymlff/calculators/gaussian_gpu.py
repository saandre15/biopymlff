import os

from biopymlff.calculators.gaussian import Gaussian

from ..cluster.slurm import get_gpu_count, get_cpu_count


class Gaussian_GPU(Gaussian):
    
    def __init__(self, label=None, xc=None, basis=None, **kwargs):
        
        super(Gaussian, self).__init__(
            label=None,
            xc=xc,
            basis=basis,
            **kwargs
        )
        
    def write_input(self, atoms, properties=None, system_changes=None):
        super().write_input(atoms, properties, system_changes)
        _path = os.path.join(self.directory, self.label + ".com")
        content = None
        with open(_path, "r") as file:
            content = file.read()

        content = "%gpucpu=" + "0-" + str(get_gpu_count() - 1) + "=" + "0-" + str(get_gpu_count() - 1) + "\n" + content # Adds GPU support
        content = "%cpu=" + "0-" + str(get_cpu_count() - 1) + "\n" + content
        vals = list(filter(lambda val: False if "TV" in val else True, content.split("\n"))) # Removes Periodic Boundary Condition
        content = ""
        for line in vals:
            content += line + "\n"
        content = content + "\n"
        
        with open(_path, "w") as file:
            file.write(content)