import os
from openmm import Platform
from openff.units import Quantity, unit
from openmm import unit as openmm_unit
from openff.toolkit import ForceField, Molecule, Topology
import numpy as np

def get_platform():
    os_platform = os.getenv('PLATFORM')
    if os_platform:
        platform = Platform.getPlatformByName(os_platform)
    else:
        # work out the fastest platform
        speed = 0
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            # print(p.getName(), p.getSpeed())
            if p.getSpeed() > speed:
                platform = p
                speed = p.getSpeed()

    print('Using platform', platform.getName())

    # if it's GPU platform set the precision to mixed
    if platform.getName() == 'CUDA' or platform.getName() == 'OpenCL':
        platform.setPropertyDefaultValue('Precision', 'mixed')
        print('Set precision for platform', platform.getName(), 'to mixed')

    return platform


def insert_molecule_and_remove_clashes(
    topology: Topology,
    insert: Molecule,
    radius: Quantity = 1.5 * unit.angstrom,
    keep: list[Molecule] = [],
) -> Topology:
    """
    Add a molecule to a copy of the topology, removing any clashing molecules.

    The molecule will be added to the end of the topology. A new topology is
    returned; the input topology will not be altered. All molecules that
    clash will be removed, and each removed molecule will be printed to stdout.
    Users are responsible for ensuring that no important molecules have been
    removed; the clash radius may be modified accordingly.

    Parameters
    ==========
    top
        The topology to insert a molecule into
    insert
        The molecule to insert
    radius
        Any atom within this distance of any atom in the insert is considered
        clashing.
    keep
        Keep copies of these molecules, even if they're clashing
    """
    # We'll collect the molecules for the output topology into a list
    new_top_mols = []
    # A molecule's positions in a topology are stored as its zeroth conformer
    insert_coordinates = insert.conformers[0][:, None, :]
    for molecule in topology.molecules:
        if any(keep_mol.is_isomorphic_with(molecule) for keep_mol in keep):
            new_top_mols.append(molecule)
            continue
        molecule_coordinates = molecule.conformers[0][None, :, :]
        diff_matrix = molecule_coordinates - insert_coordinates

        # np.linalg.norm doesn't work on Pint quantities 😢
        working_unit = unit.nanometer
        distance_matrix = (
            np.linalg.norm(diff_matrix.m_as(working_unit), axis=-1) * working_unit
        )

        if distance_matrix.min() > radius:
            # This molecule is not clashing, so add it to the topology
            new_top_mols.append(molecule)
        else:
            print(f"Removed {molecule.to_smiles()} molecule")

    # Insert the ligand at the end
    new_top_mols.append(insert)

    # This pattern of assembling a topology from a list of molecules
    # ends up being much more efficient than adding each molecule
    # to a new topology one at a time
    new_top = Topology.from_molecules(new_top_mols)

    # Don't forget the box vectors!
    new_top.box_vectors = topology.box_vectors
    return new_top



