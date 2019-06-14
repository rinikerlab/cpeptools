

from simtk.unit import *
from simtk.openmm.app import *
from simtk.openmm import *
def simulate_vacuum(mol,output_name, num_frames = 10, which_conf = 0, epsilon = 4):
    """
        writes out every 10 ps
    """

    import mlddec
    models  = mlddec.load_models(epsilon)
    charges = mlddec.get_charges(mol, models)


    from openforcefield.utils.toolkits import RDKitToolkitWrapper, ToolkitRegistry
    from openforcefield.topology import Molecule, Topology
    from openforcefield.typing.engines.smirnoff import ForceField
    # from openforcefield.typing.engines.smirnoff.forcefield import PME

    import parmed
    import numpy as np

    forcefield = ForceField('smirnoff99Frosst.offxml')


    molecule = Molecule.from_rdkit(mol)
    molecule.partial_charges = Quantity(np.array(charges), elementary_charge)
    topology = Topology.from_molecules(molecule)
    openmm_system = forcefield.create_openmm_system(topology, charge_from_molecules= [molecule])


    conf = mol.GetConformer(which_conf)
    positions = Quantity(np.array([np.array(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]), angstroms)

    structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), openmm_system, positions)




    system = structure.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)

    integrator = LangevinIntegrator(1*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(structure.topology, system, integrator)
    simulation.context.setPositions(structure.positions)


    simulation.minimizeEnergy(maxIterations = 70)

    # step = 5000
    step = 5
    simulation.reporters.append(PDBReporter(output_name, step))
    simulation.step(step * num_frames)
