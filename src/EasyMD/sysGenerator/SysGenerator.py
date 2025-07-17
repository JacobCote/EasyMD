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
from EasyMD.utils.utils import _formatIndex, writeFooter, PDBwrite_all,deletePcap

import mdtraj


class SysGenerator:
    """
    System generator for molecular dynamics simulations.
    
    This class handles the preparation of molecular systems for MD simulations, including
    protein structure fixing, ligand parameterization, solvation, and force field assignment.
    It supports both new simulations and restart from previous states, with options for
    explicit solvation or implicit solvent (GBIS) models.
    
    Attributes:
        config: Configuration object containing simulation parameters and settings.
        forcefield_kwargs (dict): Force field parameters including constraints and hydrogen mass.
        last_state (int): Current state counter for trajectory numbering.
        modeller: OpenMM Modeller object with final system topology and positions.
        system: OpenMM System object with forces and parameters ready for simulation.
    
    Example:
        >>> sys_gen = SysGenerator(config)
        >>> modeller = sys_gen.modeller  # Access prepared system
        >>> system = sys_gen.system     # Access force field system
        
    Note:
        The class automatically determines whether to prepare a new system or restart
        from a previous state based on the config.restart parameter.
    """
    def __init__(self, config):
        """
        Initialize the SysGenerator and prepare the molecular system.
        
        This constructor automatically determines whether to create a new system or
        restart from a previous state based on the configuration. It sets up force field
        parameters and prepares the complete system ready for simulation.
        
        Args:
            config: Configuration object containing all simulation parameters including
                   protein/ligand files, force fields, solvation settings, and output options.
                   
        Attributes Set:
            config: Stores the configuration object
            forcefield_kwargs (dict): Standard force field parameters with HBond constraints,
                                    rigid water, no CM motion removal, and 4 amu hydrogen mass
            last_state (int): Initialized to 0 for new simulations
            modeller: OpenMM Modeller object with prepared system
            system: OpenMM System object with force field parameters
        """
        print("SYSGENERATOR -----------------------------------")
        self.config = config
        self.forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu}
        self.last_state = 0
        
        if config.restart != None:
            print("RESTART")
            modeller, system = self._restartSetup()
        else:
            print("NO RESTART")
            modeller, system = self._setup()
        self.modeller = modeller
        self.system = system


    def _setup(self):
        """
        Set up a new molecular dynamics system from scratch.
        
        This method handles the complete preparation of a new MD system including:
        - Output directory creation and management
        - System preparation (protein-only or protein-ligand complex)
        - Force field assignment and parameterization
        
        Returns:
            Tuple[Modeller, System]: OpenMM Modeller and System objects ready for simulation.
            
        Side Effects:
            - Creates output directory (auto-numbered if not specified)
            - Calls either _prep_prot() or _prep_complex() based on ligand presence
        """
        if self.config.outdir is None:

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
            out_dir = self.config.outdir.rstrip('/')
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


    def _restartSetup(self):
        """
        Set up system for restart simulation from previously saved state.
        
        This method loads the restart configuration and prepares the system based on
        the original simulation parameters. It automatically detects whether the
        original simulation included a ligand by checking for ligand.sdf file.
        
        Returns:
            Tuple[Modeller, System]: OpenMM Modeller and System objects ready for restart.
            
        Raises:
            FileNotFoundError: If restart_setup.yml is not found in restart directory.
            KeyError: If required parameters are missing from restart setup file.
            
        Side Effects:
            - Loads restart_setup.yml configuration file
            - Calls either _prep_restart() or _prep_restart_ligand() based on ligand presence
        """
        out_dir = self.config.restart

        with open(out_dir+'/'+'restart_setup.yml') as f:
            setup = yaml.load(f, Loader=yaml.FullLoader)
            

        print('Restarting from state files in directory', self.config.restart, 'with setup file : ' )
    
        #check if ligand.sdf is present in the restart directory 
        if not os.path.isfile(out_dir+'/'+'ligand.sdf'):
            print('No ligand.sdf file found in the restart directory, preparing system without ligand ...')
            modeller, system = self._prep_restart(setup=setup, outdir=out_dir, forcefield_kwargs=self.forcefield_kwargs)
            
        else:
            print('Ligand.sdf file found in the restart directory, preparing system with ligand ...')
            modeller, system = self._prep_restart_ligand(setup=setup, outdir=out_dir, forcefield_kwargs=self.forcefield_kwargs)
        return modeller, system

        

    def _prep_restart_ligand(self, setup: dict, outdir: str, forcefield_kwargs: dict):
        """
        Prepare system for restart simulation with ligand present.
        
        This method loads the restart model PDB file and ligand SDF file, then creates
        a new system with the appropriate force fields for protein-ligand simulations.
        It handles both solvated and implicit solvent (GBIS) systems based on the
        original simulation setup.
        
        Args:
            setup (dict): Setup parameters from the original simulation including
                         force field specifications and solvation settings.
            outdir (str): Output directory containing restart files (restart_model.pdb, ligand.sdf).
            forcefield_kwargs (dict): Force field parameters for system generation.
        
        Returns:
            Tuple[Modeller, System]: OpenMM Modeller and System objects ready for restart.
        
        Raises:
            FileNotFoundError: If restart_model.pdb or ligand.sdf files are not found.
            KeyError: If required force field parameters are missing from setup.
        
        Note:
            For solvated systems, uses explicit water force fields.
            For implicit systems, uses GBIS with NoCutoff nonbonded method.
        """
        # load pdb file
        pdb = PDBFile(outdir+'/'+'restart_model.pdb')
        ligand_mol = Molecule.from_file(outdir+'/'+'ligand.sdf')
        protein_force_field = setup['protein_force_field']
        water_force_field = setup['water_force_field']
        ligand_force_field = setup['ligand_force_field']
        modeller = Modeller(pdb.topology, pdb.positions)
    

        
        if setup['solvate']:
            print('generating system with solvent...')
        
            system_generator = SystemGenerator(
            forcefields=[protein_force_field, water_force_field],
            small_molecule_forcefield=ligand_force_field,
            molecules=[ligand_mol],
            forcefield_kwargs=forcefield_kwargs)
            
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
            forcefield_kwargs=forcefield_kwargs,
            nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
        
        return modeller, system




        
    def _prep_restart(self, setup: dict, outdir: str, forcefield_kwargs: dict):
        """
        Prepare system for restart simulation without ligand (protein-only).
        
        This method loads the restart model PDB file and creates a new system with
        the appropriate force fields for protein-only simulations. It handles both
        solvated and implicit solvent (GBIS) systems based on the original simulation setup.
        
        Args:
            setup (dict): Setup parameters from the original simulation including
                         force field specifications and solvation settings.
            outdir (str): Output directory containing restart files (restart_model.pdb).
            forcefield_kwargs (dict): Force field parameters for system generation.
        
        Returns:
            Tuple[Modeller, System]: OpenMM Modeller and System objects ready for restart.
        
        Raises:
            FileNotFoundError: If restart_model.pdb file is not found.
            KeyError: If required force field parameters are missing from setup.
        
        Note:
            For solvated systems, uses explicit water force fields.
            For implicit systems, uses GBIS with NoCutoff nonbonded method.
        """
        # load pdb file
        pdb = PDBFile(outdir+'/'+'restart_model.pdb')
        protein_force_field = setup['protein_force_field']
        water_force_field = setup['water_force_field']
        modeller = Modeller(pdb.topology, pdb.positions)
    

        
        if setup['solvate']:
            print('generating system with solvent...')
        
            system_generator = SystemGenerator(
            forcefields=[protein_force_field, water_force_field],
            forcefield_kwargs=forcefield_kwargs)

            system = system_generator.create_system(modeller.topology)
        
        
        else :
            
            pdb = mdtraj.load(outdir+'/'+'restart_model.pdb')
            topology = pdb.topology.to_openmm()
            print('generating system without solvent...')
            system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
        
            forcefield_kwargs=forcefield_kwargs,
            nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            
            system = system_generator.create_system(topology)
        
        return modeller, system


        




    def _prep_prot(self, pdb_in, list_of_molecules_to_remove,
                    solvate, protein_force_field, water_force_field,
                    water_model, positive_ion,
                    negative_ion, ionic_strength, no_neutralize, padding, ph,
                    outdir, forcefield_kwargs):
        """
        Prepare protein-only system for molecular dynamics simulation.
        
        This method performs comprehensive protein preparation including structure fixing,
        molecule removal, solvation (if requested), and force field assignment. It handles
        both explicit solvent and implicit solvent (GBIS) systems.
        
        Args:
            pdb_in (str): Path to input PDB file containing the protein structure.
            list_of_molecules_to_remove (list): List of molecule names to remove from structure.
            solvate (bool): Whether to add explicit solvent box.
            protein_force_field (str): Protein force field specification (e.g., 'amber14-all.xml').
            water_force_field (str): Water force field specification (e.g., 'amber/tip3p_standard.xml').
            water_model (str): Water model to use for solvation (e.g., 'tip3p').
            positive_ion (str): Positive ion type for neutralization (e.g., 'Na+').
            negative_ion (str): Negative ion type for neutralization (e.g., 'Cl-').
            ionic_strength (float): Ionic strength for solvation (M).
            no_neutralize (bool): Whether to skip system neutralization.
            padding (float): Solvent box padding around protein (Å).
            ph (float): pH for protonation state assignment.
            outdir (str): Output directory for generated files.
            forcefield_kwargs (dict): Force field parameters and constraints.
        
        Returns:
            Tuple[Modeller, System]: OpenMM Modeller and System objects ready for simulation.
        
        Side Effects:
            - Creates prot_receptor.pdb: Cleaned protein structure
            - Creates solvated_complex.pdb: Solvated system (if solvate=True)
            - Creates restart_model.pdb: Final system ready for simulation
            
        Raises:
            FileNotFoundError: If input PDB file is not found.
            Exception: If protein fixing or solvation fails.
            
        Note:
            The method automatically adds missing residues, atoms, and hydrogens,
            removes specified molecules (including water), and optionally solvates
            the system with specified parameters.
        """
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



                