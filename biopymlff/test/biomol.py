from general import General_Test

class Biomol_Test(General_Test):
    
    @classmethod
    def setUpClass(cls):
        cls.source_traj = []
        cls.target_traj = []

        source_traj = []
        target_traj = []

        mol.calc = self.source
        dynamics = Langevin(atoms=mol, timestep=0.0002, temperature_K=500)
        collect_data = lambda: source_traj.append(ase_atoms_to_pytraj_atoms(mol.copy()))
        dynamics.attach(collect_data, interval=1)
        dynamics.run(steps=100)

        mol = read_proteindatabank(testing_params["pdb_file"])
        mol.calc = self.target
        dynamics = Langevin(atoms=mol, timestep=0.0002, temperature_K=500)
        collect_data = lambda: target_traj.append(ase_atoms_to_pytraj_atoms(mol.copy()))
        dynamics.attach(collect_data, interval=1)
        dynamics.run(steps=100)

        super().setUpClass(cls)

    def get_source_frames(dt: float, frame_count: int):
        # Perform ASE->PyTraj Atom Conversion
        pass

    def get_target_frames(dt: float, frame_count: int):
        pass