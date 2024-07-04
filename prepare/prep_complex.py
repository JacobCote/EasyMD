from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
from openmm import app, unit
from openmm.app import PDBFile, Modeller
from openff.toolkit import  Molecule
from pdbfixer import PDBFixer
from utils.openbabel_charge import get_charges
from utils.utils import _formatIndex, writeFooter, PDBwrite_all,deletePcap
import mdtraj


def prep_complex(pdb_in,list_of_molecules_to_remove,lig_name,
                 solvate,protein_force_field,water_force_field,
                 ligand_force_field,water_model,positive_ion,
                 negative_ion,ionic_strength,no_neutralize,padding,ph,
                 outdir,forcefield_kwargs):
    
    # load and fix the protein
        
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    print('Residues:', fixer.missingResidues)
    print('Atoms:', fixer.missingAtoms)
    print('Terminals:', fixer.missingTerminals)
    print('Non-standard:', fixer.nonstandardResidues)
    print('Changing non-standard residues to standard residues...')
    fixer.replaceNonstandardResidues()
    
    

    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    
    
    # remove molecules from pdb and save it in a new file
    list_of_molecules_to_remove += ['HOH','WAT',lig_name]
    modeller = Modeller(fixer.topology, fixer.positions)
    print('Before remove water, System has %d atoms' % modeller.topology.getNumAtoms())
    for mol in list_of_molecules_to_remove:
       
        toDelete = []
        for res in modeller.topology.residues():
            if res.name == mol:
                toDelete.append(res)
                print('Deleting', res)
        modeller.delete(toDelete)
    print('After remove water, System has %d atoms' % modeller.topology.getNumAtoms())
    
    PDBwrite_all(modeller, outdir+'/'+'prot_receptor.pdb')


    print('Done')
    # extract ligan from pdb and save it in a new file
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    print('Residues:', fixer.missingResidues)
    print('Atoms:', fixer.missingAtoms)
    print('Terminals:', fixer.missingTerminals)
    print('Non-standard:', fixer.nonstandardResidues)

    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    modeller = Modeller(fixer.topology, fixer.positions)
    toDelete = []
    for res in modeller.topology.residues():
        if res.name != lig_name:
            toDelete.append(res)
    modeller.delete(toDelete)

    with open(outdir+'/'+'ligand.pdb', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)
    # calculate charges with openbabel anf export as sdf file
    get_charges(outdir+'/'+'ligand.pdb',outdir+'/'+'ligand.sdf','pdb','sdf')

        

        
    # calculate charges with openbabel adn export as sdf file
    
    
    # load the ligand and the protein
    protein_pdb = PDBFile(outdir+'/'+'prot_receptor.pdb')
    ligand_mol = Molecule.from_file(outdir+'/'+'ligand.sdf')


    # The topology is described in the openforcefield API
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
    print('System has %d atoms' % modeller.topology.getNumAtoms())
    
    # The topology is described in the openforcefield API
    print('Adding ligand...')
    lig_top = ligand_mol.to_topology()
    modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
    print('System has %d atoms' % modeller.topology.getNumAtoms())
    # write to pdb
    with open(outdir+'/'+'complex.pdb', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    
    # replace UNK with the ligand name in the pdb file 
    with open(outdir+'/'+'complex.pdb', 'r') as file :
        filedata = file.read() 
        filedata = filedata.replace('UNK',lig_name)
        
    PDBwrite_all(modeller, outdir+'/'+'complex.pdb')

        
    
    
    if solvate:
        print('generating system with solvent...')
        modeller = deletePcap(modeller)
       
        system_generator = SystemGenerator(
        forcefields=[protein_force_field, water_force_field],
        small_molecule_forcefield=ligand_force_field,
        molecules=[ligand_mol],
        forcefield_kwargs=forcefield_kwargs)
        
        print('Adding solvent...')
        # we use the 'padding' option to define the periodic box.
        # we just create a box that has a 10A (default) padding around the complex.
        modeller.addSolvent(system_generator.forcefield, model=water_model, padding=padding * unit.angstroms,
                        positiveIon=positive_ion, negativeIon=negative_ion,
                        ionicStrength=ionic_strength * unit.molar, neutralize=not no_neutralize)
        print('System has %d atoms' % modeller.topology.getNumAtoms())

        with open(outdir+'/'+'solvated_complex.pdb', 'w') as outfile:
            PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
        print('IMPORTANT: The solvated complex is saved as solvated_complex.pdb in the output directory, please check it before proceeding')
        pdb = PDBFile(outdir+'/'+'solvated_complex.pdb')
        modeller = Modeller(pdb.topology, pdb.positions)
        system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
    
    
    else :
        #loading the pdb pdb file with mdtraj not forcing me to have a periodic system ???
        pdb = mdtraj.load(outdir+'/'+'complex.pdb')
        modeller = Modeller(pdb.topology.to_openmm(), pdb.xyz[0])
        # remove pcap from adn
        modeller = deletePcap(modeller)
        
        print('generating system without solvent...')
        system_generator = SystemGenerator(
        forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
        small_molecule_forcefield=ligand_force_field,
        molecules=[ligand_mol],
        forcefield_kwargs=forcefield_kwargs,
        nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
        )
        system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
        
    # write to pdb modeller as restart_model.pdb
    
    PDBwrite_all(modeller, outdir+'/'+'restart_model.pdb')
   
    return modeller, system
