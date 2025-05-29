from openmm import app, unit
import yaml
import os
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator

from openmm.app import PDBFile,  Modeller
from openff.toolkit import Molecule
import yaml
import mdtraj

from pdbfixer import PDBFixer
from utils.utils import _formatIndex, writeFooter, PDBwrite_all,deletePcap
import mdtraj



class SysGenerator():
    def __init__(self,config):
        self.config = config
        self.forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
        self.last_state = 0
        if config.restart :
            self._restartSetup()
        else :
            self._setup()



    def _setup(self,):
        if self.config.output is None:

            exist = True
            i = 0
            while exist:
                out_dir = 'out_'+str(i)
                exist = os.path.isdir(out_dir)
                if not exist:
                    os.mkdir(out_dir)
                i += 1

        else:
            # strip the trailing slash
            out_dir = self.config.output.rstrip('/')
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)
            



        pdb_in = self.config.protein
        mol_in = self.config.ligand
        output_traj_dcd = f'output_traj_{self.last_state}.dcd'
        num_steps = self.config.steps
        reporting_interval = self.config.interval
        temperature = self.config.temperature * unit.kelvin
        equilibration_steps = self.config.equilibration_steps

        if self.config.ligand is None:
            modeller, system = self._prep_prot(
                pdb_in = pdb_in,
                list_of_molecules_to_remove=self.config.remove,
                solvate=self.config.solvate,
                protein_force_field=self.config.protein_force_field,
                water_force_field=self.config.water_force_field,
                water_model=self.config.water_model,
                positive_ion=self.config.positive_ion,
                negative_ion=self.config.negative_ion,
                ionic_strength=self.config.ionic_strength,
                no_neutralize=self.config.no_neutralize,
                padding=self.config.padding,
                ph=self.config.ph,
                outdir=out_dir,
                forcefield_kwargs=self.forcefield_kwargs
            )
            
        else:
            modeller,system = self._prep_complex(
                pdb_in = pdb_in,
                list_of_molecules_to_remove=self.config.remove,
                lig_name=self.config.ligand,
                solvate=self.config.solvate,
                protein_force_field=self.config.protein_force_field,
                water_force_field=self.config.water_force_field,
                ligand_force_field=self.config.ligand_force_field,
                water_model=self.config.water_model,
                positive_ion=self.config.positive_ion,
                negative_ion=self.config.negative_ion,
                ionic_strength=self.config.ionic_strength,
                no_neutralize=self.config.no_neutralize,
                padding=self.config.padding,
                ph=self.config.ph,
                outdir=out_dir,
                forcefield_kwargs=self.forcefield_kwargs
            )
        return modeller, system


    def _restartSetup(self,):
        out_dir = self.restart_dir
        with open(out_dir+'/'+'restart_setup.yml') as f:
            setup = yaml.load(f, Loader=yaml.FullLoader)

        print('Restarting from state files in directory', self.config.restart_dir, 'with setup file : ' )
    
        #check if ligand.sdf is present in the restart directory 
        if not os.path.isfile(out_dir+'/'+'ligand.sdf'):
            print('No ligand.sdf file found in the restart directory, preparing system without ligand ...')
            modeller, system = self._prep_restart(setup=setup, outdir=out_dir, forcefield_kwargs=self.forcefield_kwargs)
            
        else:
            print('Ligand.sdf file found in the restart directory, preparing system with ligand ...')
            modeller, system = self._prep_restart_ligand(setup=setup, outdir=out_dir, forcefield_kwargs=self.forcefield_kwargs)
        return modeller, system

        

    def _prep_restart_ligand(self,):
        '''
        prepare system for restart simulation with ligand
        :param setup: dict, setup parameters for the simulation wich is created when the simulation is started for the first time
        :param outdir: str, output directory where the restart files are stored
        :param forcefield_kwargs: dict, forcefield parameters
        :return: Tuple[Modeller, SystemGenerator], Modeller object and SystemGenerator object
        '''
        # load pdb file
        pdb = PDBFile(self.config.outdir+'/'+'restart_model.pdb')
        ligand_mol = Molecule.from_file(self.config.outdir+'/'+'ligand.sdf')
        protein_force_field = self.config.protein_force_field
        water_force_field = self.config.water_force_field
        ligand_force_field = self.config.ligand_force_field
        modeller = Modeller(pdb.topology, pdb.positions)
    

        
        if self.config.solvate:
            print('generating system with solvent...')
        
            system_generator = SystemGenerator(
            forcefields=[protein_force_field, water_force_field],
            small_molecule_forcefield=ligand_force_field,
            molecules=[ligand_mol],
            forcefield_kwargs=self.forcefield_kwargs)
            
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
        
        
        else :
            
            '''
            pdb = mdtraj.load(outdir+'/'+'complex.pdb')
            topology = pdb.topology.to_openmm()
            '''
            print('generating system without solvent...')
            system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
            small_molecule_forcefield=ligand_force_field,
            molecules=[ligand_mol],
            forcefield_kwargs=self.forcefield_kwargs,
            nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
        
        return modeller, system




    def _prep_restart(self ):
        '''
        prepare system for restart simulation without ligand
        :param setup: dict, setup parameters for the simulation wich is created when the simulation is started for the first time
        :param outdir: str, output directory where the restart files are stored
        :param forcefield_kwargs: dict, forcefield parameters
        :return: Tuple[Modeller, SystemGenerator], Modeller object and SystemGenerator object
        '''
        # load pdb file
        pdb = PDBFile(self.config.outdir+'/'+'restart_model.pdb')
        protein_force_field = self.config.protein_force_field
        water_force_field = self.config.water_force_field
        modeller = Modeller(pdb.topology, pdb.positions)
    

        
        if self.config.solvate:
            print('generating system with solvent...')
        
            system_generator = SystemGenerator(
            forcefields=[protein_force_field, water_force_field],
            forcefield_kwargs=self.forcefield_kwargs)

            system = system_generator.create_system(modeller.topology)
        
        
        else :
            
            pdb = mdtraj.load(self.config.outdir+'/'+'restart_model.pdb')
            topology = pdb.topology.to_openmm()
            print('generating system without solvent...')
            system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
        
            forcefield_kwargs=self.forcefield_kwargs,
            nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            
            system = system_generator.create_system(topology)
        
        return modeller, system
    


    def _prep_prot(self,pdb_in,list_of_molecules_to_remove,
                    solvate,protein_force_field,water_force_field,
                    water_model,positive_ion,
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
        list_of_molecules_to_remove += ['HOH','WAT']
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

        # The topology is described in the openforcefield API
        #modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
        print('System has %d atoms' % modeller.topology.getNumAtoms())
        
        
        if solvate:
            modeller = deletePcap(modeller)
            print('generating system with solvent...')
        
            system_generator = SystemGenerator(
            forcefields=[protein_force_field, water_force_field],
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
            system = system_generator.create_system(modeller.topology)
        
        
        else :
            # loading the pdb pdb file with mdtraj not forcing me to have a periodic system ???
            pdb = mdtraj.load(outdir+'/'+'prot_receptor.pdb')
            modeller = Modeller(pdb.topology.to_openmm(), pdb.xyz[0])
            modeller = deletePcap(modeller)
            topology = modeller.topology
            
            #topology = pdb.topology.to_openmm()
            print('generating system without solvent...')
            system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
            forcefield_kwargs=forcefield_kwargs,
            nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            system = system_generator.create_system(topology)
            print('System crated without solvent')
            
        # write to pdb modeller as restart_model.pdb
    
        PDBwrite_all(modeller, outdir+'/'+'restart_model.pdb')
        
        
        return modeller, system



                