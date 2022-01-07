import os
import uuid

from ase.calculators.gaussian import Gaussian as G

from biopymlff.util.getenv import getenv

class Gaussian(G):
    
    def __init__(self, label=None, xc=None, basis=None, **kwargs):
        general_params = getenv()['general']
        gaussian_params = getenv()['gaussian']

        directory = general_params['scratch_dir']
        if "$" in directory: directory = os.environ[directory.replace("$", "")]

        super(G, self).__init__(
            label=label,
            directory=os.path.join(directory, "gaussian", label + "_" + uuid.uuid1().hex),
            # mem=gaussian_params['memory'] \
            #     if gaussian_params['memory'] != "auto" \
            #     else str((os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024.**3)) - 5) + "GB",
            mem="10GB",
            chk=gaussian_params['checkpoint_file'] if gaussian_params['checkpoint_file'] else label + ".chk",
            save=None,
            xc=xc,
            basis=basis,
            # scf='noincfock,novaracc,fermi,maxcycle=3000,ndamp=64,xqc',
            # int='acc2e=12',
            **kwargs
        )

    