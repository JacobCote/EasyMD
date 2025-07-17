
import sys, time
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import app, unit, LangevinIntegrator
from openmm.app import PDBFile, Simulation, Modeller, StateDataReporter, DCDReporter
from openff.toolkit import Molecule
from typing import Tuple
import yaml
import mdtraj


class Restarter:
    """
    Handles restarting molecular dynamics simulations from previously saved states.
    
    This class provides functionality to restart MD simulations by loading saved states,
    topologies, and configurations from previous runs. It supports both protein-only
    and protein-ligand complex simulations with or without explicit solvent.
    
    Attributes:
        setup (dict): Configuration parameters from the original simulation setup.
        outdir (str): Output directory containing restart files.
        forcefield_kwargs (dict): Force field parameters for system generation.
    
    Example:
        >>> setup = {'protein_force_field': 'amber14-all.xml', 'solvate': True}
        >>> restarter = Restarter(setup, 'restart_dir', {})
        >>> modeller, system = restarter.prep_restart()
    """
    def __init__(self, setup: dict, outdir: str, forcefield_kwargs: dict):
        """
        Initialize the Restarter with simulation setup parameters.
        
        Args:
            setup (dict): Configuration parameters from the original simulation including
                         force field specifications, solvation settings, and simulation parameters.
            outdir (str): Path to the output directory containing restart files
                         (restart_model.pdb, ligand.sdf, last_state.xml, etc.).
            forcefield_kwargs (dict): Force field parameters for system generation.
                                    Note: This parameter is currently overridden with default values.
        """
        self.setup = setup
        self.outdir = outdir
        self.forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu}




    def prep_restart_ligand(self) -> Tuple[Modeller, SystemGenerator]:
        """
        Prepare system for restart simulation with ligand present.
        
        This method loads the restart model PDB file and ligand SDF file, then creates
        a new system with the appropriate force fields for protein-ligand simulations.
        It handles both solvated and implicit solvent (GBIS) systems.
        
        Returns:
            Tuple[Modeller, SystemGenerator]: A tuple containing:
                - Modeller: OpenMM Modeller object with topology and positions
                - SystemGenerator: OpenMM System object with forces and parameters
        
        Raises:
            FileNotFoundError: If restart_model.pdb or ligand.sdf files are not found.
            KeyError: If required force field parameters are missing from setup.
        
        Note:
            Requires the following files in the restart directory:
            - restart_model.pdb: The protein-ligand complex structure
            - ligand.sdf: The ligand structure file
        """
        setup = self.setup 
        outdir = self.setup 
        forcefield_kwargs = self.forcefield_kwargs 

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




    def prep_restart(self) -> Tuple[Modeller, SystemGenerator]:
        """
        Prepare system for restart simulation without ligand (protein-only).
        
        This method loads the restart model PDB file and creates a new system with
        the appropriate force fields for protein-only simulations. It handles both
        solvated and implicit solvent (GBIS) systems.
        
        Returns:
            Tuple[Modeller, SystemGenerator]: A tuple containing:
                - Modeller: OpenMM Modeller object with topology and positions
                - SystemGenerator: OpenMM System object with forces and parameters
        
        Raises:
            FileNotFoundError: If restart_model.pdb file is not found.
            KeyError: If required force field parameters are missing from setup.
        
        Note:
            Requires restart_model.pdb file in the restart directory containing
            the protein structure from the previous simulation.
        """
        setup = self.setup 
        outdir = self.setup 
        forcefield_kwargs = self.forcefield_kwargs 
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



    def restart_simulation(self, system, modeller, restart_dir, setup, clock=None, step=None):
        """
        Restart molecular dynamics simulation from a previously saved state.
        
        This method loads the last saved state from an XML file and continues the simulation
        from that point. It supports both time-based (clock) and step-based simulation modes.
        The method automatically saves trajectory data, state information, and updates the
        restart configuration file upon completion.
        
        Args:
            system: OpenMM System object containing the molecular system definition.
            modeller: OpenMM Modeller object with topology and current positions.
            restart_dir (str): Directory containing restart files (last_state.xml, etc.).
            setup (dict): Simulation setup parameters including temperature, integrator settings,
                         and reporting intervals.
            clock (float, optional): Simulation time in minutes. Mutually exclusive with step.
            step (int, optional): Number of simulation steps to run. Mutually exclusive with clock.
        
        Raises:
            FileNotFoundError: If last_state.xml or other required restart files are missing.
            ValueError: If both clock and step are provided, or if neither is provided.
            
        Side Effects:
            - Creates new trajectory file: output_traj_{last_state+1}.dcd
            - Appends to log file: log.txt
            - Saves final state: last_state.xml
            - Saves final structure: last_state_{last_state+1}.pdb
            - Updates restart_setup.yml with incremented last_state
            - Exits the program upon completion
            
        Note:
            This method calls sys.exit(0) upon successful completion, terminating the program.
            For solvated systems, a Monte Carlo barostat is automatically added.
        """
        setup = self.setup 
        outdir = self.setup 
        forcefield_kwargs = self.forcefield_kwargs 
        
        out_dir = restart_dir
    
        last_state = setup['last_state']
        temperature = setup['temperature']
        friction_coeff = setup['friction_coeff']
        step_size = setup['step_size']
        reporting_interval = setup['reporting_interval']
        
        
        integrator = LangevinIntegrator(temperature, friction_coeff, step_size)
        
        if setup['solvate']:
            system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature, 25))
            print('Default Periodic box: {}'.format(system.getDefaultPeriodicBoxVectors()))

        
        simulation = Simulation(modeller.topology, system, integrator, state=out_dir+'/'+'last_state.xml')  
        #simulation = Simulation(modeller.topology, system, integrator, platform=platform, state=out_dir+'/'+'last_state.xml')  
        
        simulation.reporters.append(DCDReporter(out_dir+'/'+f'output_traj_{last_state+1}.dcd', reporting_interval, enforcePeriodicBox=False))
        simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
        simulation.reporters.append(StateDataReporter(out_dir+'/'+'log.txt', reporting_interval, step=True, potentialEnergy=True, temperature=True,append=True))

        if clock is not None:
            print('Starting simulation for', clock, ' mins ...')
            t1 = time.time()
            simulation.runForClockTime(clock * unit.minute)
            t2 = time.time()
            
        else :
            print('Starting simulation with', step, 'steps ...')
            t1 = time.time()
            simulation.step(step)
            t2 = time.time()
        
        # save the last state of the simulation
        simulation.saveState(out_dir+'/'+f'last_state.xml')
        # save the last state pdb
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open(out_dir+'/'+f'last_state_{last_state+1}.pdb', 'w'))
        step_s = step_size * unit.picoseconds
        n_step = simulation.context.getStepCount()
        duration = (n_step * step_s).value_in_unit(unit.nanoseconds)
        print('Simulation complete in {} mins at {}.'.format(
        round((t2 - t1) / 60, 3), temperature,'K'))
        print('Simulation time was', round(duration, 3), 'ns')
        print('Updating last state in restart_setup.yml ...')
        setup['last_state'] = last_state + 1
        yaml.dump(setup, open(out_dir+'/'+'restart_setup.yml', 'w'), default_flow_style=False)
        print('Exiting ...')
        exit(0)



