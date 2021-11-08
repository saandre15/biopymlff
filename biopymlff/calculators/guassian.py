from ase.calculators.gaussian import Gaussian as G

class Gaussian(G):
    
    command = 'module use /work2/01114/jfonner/frontera/modulefiles; module load gaussian; ' + G.command