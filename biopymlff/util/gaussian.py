import os

from getenv import getenv

from ase.calculators.gaussian import Gaussian


def get_gaussian(label: str, xc=None, basis=None):
    general_params = getenv()['general']
    gaussian_params = getenv()['gaussian']
    return Gaussian(
        label=label,
        directory=general_params['scratch_dir'] + "/gaussian/" + self.label + "_" + uuid.uuid1().hex,
        mem=gaussian_params['memory'] \
            if gaussian_params['memory'] != "auto" \
            else str((os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024.**3)) - 5) + "GB",
        chk=gaussian_params['checkpoint_file'] if gaussian_params['checkpoint_file'] else label + ".chk",
        save=None,
        xc=method,
        basis=basis,
        scf='(noincfock,novaracc,fermi,maxcycle=3000,ndamp=64,xqc)',
        int='acc2e=12'
    )