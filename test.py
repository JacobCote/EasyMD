from openmm.app import PDBFile, Modeller
import pandas as pd
import xml.etree.ElementTree as ET
from pdbfixer import PDBFixer

ph = 7.0
fixer = PDBFixer(filename='p33_adn_og.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.findNonstandardResidues()
print('Residues:', fixer.missingResidues)
print('Atoms:', fixer.missingAtoms)
print('Terminals:', fixer.missingTerminals)
print('Non-standard:', fixer.nonstandardResidues)
print('Changing non-standard residues to standard residues...')
fixer.replaceNonstandardResidues()
modeller = Modeller(fixer.topology, fixer.positions)

PDBFile.writeHeader(modeller.topology, open('test.pdb', 'w'))
PDBFile.writeFile(modeller.topology, modeller.positions, file=open('test.pdb', 'a'), keepIds=True)
PDBFile.writeFooter(modeller.topology, open('test2.pdb', 'w'))



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