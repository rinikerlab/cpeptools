
def minimise_energy_all_confs(mol, models = None, epsilon = 4, allow_undefined_stereo = True, **kwargs ):
    from simtk.unit import *
    from simtk.openmm.app import *
    from simtk.openmm import *
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    import mlddec
    import copy
    import tqdm
    mol = Chem.AddHs(mol, addCoords = True)

    if models is None:
        models  = mlddec.load_models(epsilon)
    charges = mlddec.get_charges(mol, models)

    from openforcefield.utils.toolkits import RDKitToolkitWrapper, ToolkitRegistry
    from openforcefield.topology import Molecule, Topology
    from openforcefield.typing.engines.smirnoff import ForceField
    # from openforcefield.typing.engines.smirnoff.forcefield import PME

    import parmed
    import numpy as np

    forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

    # molecule = Molecule.from_rdkit(mol, allow_undefined_stereo = True)
    molecule = Molecule.from_rdkit(mol, allow_undefined_stereo = allow_undefined_stereo)
    molecule.partial_charges = Quantity(np.array(charges), elementary_charge)
    topology = Topology.from_molecules(molecule)
    openmm_system = forcefield.create_openmm_system(topology, charge_from_molecules= [molecule])


    # conf = mol.GetConformer(0)
    # positions = Quantity(np.array([np.array(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]), angstroms)

    structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), openmm_system)


    system = structure.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)

    integrator = LangevinIntegrator(273*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(structure.topology, system, integrator)

    out_mol = copy.deepcopy(mol)
    for i in tqdm.tqdm(range(out_mol.GetNumConformers())):
        conf = mol.GetConformer(i)
        structure.coordinates =  Quantity(np.array([np.array(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]), angstroms)

        simulation.context.setPositions(structure.positions)

        simulation.minimizeEnergy()
        simulation.step(1)

        coords = simulation.context.getState(getPositions = True).getPositions(asNumpy = True).value_in_unit(angstrom)
        conf = out_mol.GetConformer(i)
        for j in range(out_mol.GetNumAtoms()):
            conf.SetAtomPosition(j, Point3D(*coords[j]))

    return out_mol

"""
def parameterise_molecule(mol, which_conf = -1, models = None, epsilon = 4, allow_undefined_stereo = True, **kwargs):
    """
    which_conf : when -1 write out multiple parmed structures, one for each conformer
    """
    from rdkit import Chem
    import mlddec
    mol = Chem.AddHs(mol, addCoords = True)

    if models is None:
        models  = mlddec.load_models(epsilon)
    charges = mlddec.get_charges(mol, models)


    from openforcefield.utils.toolkits import RDKitToolkitWrapper, ToolkitRegistry
    from openforcefield.topology import Molecule, Topology
    from openforcefield.typing.engines.smirnoff import ForceField
    # from openforcefield.typing.engines.smirnoff.forcefield import PME

    import parmed
    import numpy as np

    forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

    # molecule = Molecule.from_rdkit(mol, allow_undefined_stereo = True)
    molecule = Molecule.from_rdkit(mol, allow_undefined_stereo = allow_undefined_stereo)
    molecule.partial_charges = Quantity(np.array(charges), elementary_charge)
    topology = Topology.from_molecules(molecule)
    openmm_system = forcefield.create_openmm_system(topology, charge_from_molecules= [molecule])


    if which_conf == -1 : #TODO better design here
        for i in range(mol.GetNumConformers()):
            conf = mol.GetConformer(which_conf)
            positions = Quantity(np.array([np.array(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]), angstroms)
            structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), openmm_system, positions)
            yield structure

    else:
        conf = mol.GetConformer(which_conf)
        positions = Quantity(np.array([np.array(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]), angstroms)

        structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), openmm_system, positions)
        yield structure

def parameterise_molecule_am1bcc(mol, which_conf = 0, allow_undefined_stereo = True, **kwargs):
    from rdkit import Chem
    mol = Chem.AddHs(mol, addCoords = True)

    from openforcefield.topology import Molecule, Topology
    from openforcefield.typing.engines.smirnoff import ForceField
    # from openforcefield.typing.engines.smirnoff.forcefield import PME
    from openforcefield.utils.toolkits import AmberToolsToolkitWrapper

    import parmed
    import numpy as np

    forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

    # molecule = Molecule.from_rdkit(mol, allow_undefined_stereo = True)
    molecule = Molecule.from_rdkit(mol, allow_undefined_stereo = allow_undefined_stereo)

    molecule.compute_partial_charges_am1bcc(toolkit_registry = AmberToolsToolkitWrapper())

    topology = Topology.from_molecules(molecule)
    openmm_system = forcefield.create_openmm_system(topology, charge_from_molecules= [molecule])


    conf = mol.GetConformer(which_conf)
    positions = Quantity(np.array([np.array(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]), angstroms)

    structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), openmm_system, positions)
    return structure


def minimise_energy(structure, output_name, **kwargs):
    # structure = parameterise_molecule_am1bcc(mol, **kwargs)

    system = structure.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)

    integrator = LangevinIntegrator(273*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(structure.topology, system, integrator)
    simulation.context.setPositions(structure.positions)

    simulation.minimizeEnergy()

    simulation.reporters.append(PDBReporter(output_name, 1))
    simulation.step(1)
    return simulation

def simulate_vacuum(structure, output_name, num_frames = 10):
    """
        writes out every 10 ps
    """
    # structure = parameterise_molecule(mol)

    system = structure.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)

    integrator = LangevinIntegrator(1*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(structure.topology, system, integrator)
    simulation.context.setPositions(structure.positions)


    simulation.minimizeEnergy(maxIterations = 50)

    step = 5000
    # step = 5
    simulation.reporters.append(PDBReporter(output_name, step))
    simulation.step(step * num_frames)
"""
