import os
from openmm import Platform
from openff.units import Quantity, unit
from openmm import unit as openmm_unit
from openff.toolkit import ForceField, Molecule, Topology
import numpy as np
from openmm.app import PDBFile, Modeller

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


def writeFooter(topology, file):
    """Write out the footer for a PDB file.

    Parameters
    ----------
    topology : Topology
        The Topology defining the molecular system being written
    file : file=stdout
        A file to write the file to
    """
    # Identify bonds that should be listed as CONECT records.

    conectBonds = []
    for atom1, atom2 in topology.bonds():
        if atom1.residue.name in {'DC','DG','DT','DA'} or atom2.residue.name in {'DC','DG','DT','DA'}:
            conectBonds.append((atom1, atom2))
        elif atom1.name == 'SG' and atom2.name == 'SG' and atom1.residue.name == 'CYS' and atom2.residue.name == 'CYS':
            conectBonds.append((atom1, atom2))
    if len(conectBonds) > 0:

        # Work out the index used in the PDB file for each atom.

        atomIndex = {}
        nextAtomIndex = 0
        prevChain = None
        for chain in topology.chains():
            for atom in chain.atoms():
                if atom.residue.chain != prevChain:
                    nextAtomIndex += 1
                    prevChain = atom.residue.chain
                atomIndex[atom] = nextAtomIndex
                nextAtomIndex += 1

        # Record which other atoms each atom is bonded to.

        atomBonds = {}
        for atom1, atom2 in conectBonds:
            index1 = atomIndex[atom1]
            index2 = atomIndex[atom2]
            if index1 not in atomBonds:
                atomBonds[index1] = []
            if index2 not in atomBonds:
                atomBonds[index2] = []
            atomBonds[index1].append(index2)
            atomBonds[index2].append(index1)

        # Write the CONECT records.

        for index1 in sorted(atomBonds):
            bonded = atomBonds[index1]
            while len(bonded) > 4:
                print("CONECT%5s%5s%5s%5s" % (_formatIndex(index1, 5), _formatIndex(bonded[0], 5), _formatIndex(bonded[1], 5), _formatIndex(bonded[2], 5)), file=file)
                del bonded[:4]
            line = "CONECT%5s" % _formatIndex(index1, 5)
            for index2 in bonded:
                line = "%s%5s" % (line, _formatIndex(index2, 5))
            print(line, file=file)
    print("END", file=file)

def _formatIndex(index, places):
    """Create a string representation of an atom or residue index.  If the value is larger than can fit
    in the available space, switch to hex.
    """
    if index < 10**places:
        format = f'%{places}d'
        return format % index
    format = f'%{places}X'
    shiftedIndex = (index - 10**places + 10*16**(places-1)) % (16**places)
    return format % shiftedIndex

def PDBwrite_all(modeller, filename):
    with open(filename, 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
        
    with open(filename, "r+") as f:
        lines = f.readlines()    
        
    with open(filename, "w") as f:
        for line in lines[:-1]:
            f.write(line)   
    with open(filename, 'a') as outfile:
        writeFooter(modeller.topology, outfile)
    
    

def deletePcap(modeller):
    toDelete = []
    chains = set()
    for res in modeller.topology.residues():
        if res.name in ['DC','DT','DA','DG']:
            chains.add(res.chain)


    for res in modeller.topology.residues():
        if res.chain in chains: 
            chains.remove(res.chain)
            for atom in res.atoms():
                if atom.name in {'OP3', 'OP1', 'OP2', 'P'}:
                    toDelete.append(atom)
                    print('Deleting PCAP : ', atom)
    modeller.delete(toDelete)
    return modeller

                        


