
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
    Molecular dynamics simulation runner for restarting simulations from saved states.
    
    This class handles the continuation of MD simulations from previously saved states,
    allowing for long simulations to be run in segments or for simulations to be
    resumed after interruption. It loads the previous state and continues the simulation
    with the same parameters and conditions.
    
    Attributes:
        setup (dict): Configuration parameters from the original simulation setup.
        config: Current configuration object with simulation parameters.
        modeller: OpenMM Modeller object with system topology and positions.
        system: OpenMM System object with forces and parameters.
        temperature (unit.Quantity): Simulation temperature in Kelvin.
        equilibration_steps (int): Number of equilibration steps (not used in restart).
        step_size (unit.Quantity): Integration time step in picoseconds.
    
    Example:
        >>> setup = yaml.load(open('restart_setup.yml'))
        >>> restarter = Restarter(setup, config, modeller, system)
        >>> restarter.run()  # Continues simulation from last saved state
        
    Note:
        Requires last_state.xml file and restart_setup.yml configuration file
        from the previous simulation run.
    """
    def __init__(self, setup, config, modeller, system):
        """
        Initialize the Restarter with simulation parameters and previous setup.
        
        Args:
            setup (dict): Configuration parameters from the original simulation including
                         last_state, temperature, friction_coeff, step_size, and reporting_interval.
            config: Current configuration object with simulation parameters (steps, clock, restart directory).
            modeller: OpenMM Modeller object with system topology and positions.
            system: OpenMM System object containing force field parameters and constraints.
        """
        
        self.setup = setup
        self.setup = config.config
        
       
        
        
        self.config = config
       
        self.modeller = modeller
        self.system = system
        self.temperature = self.config.temperature * unit.kelvin
        self.equilibration_steps = self.config.equilibration_steps
        self.step_size = self.config.step_size * unit.picoseconds
      
        
    
    def run(self):
        """
        Execute the restart molecular dynamics simulation from a previously saved state.
        
        This method loads the last saved state from an XML file and continues the simulation
        from that exact point, maintaining continuity with the previous run. It supports both
        time-based (clock) and step-based simulation modes, and automatically handles
        trajectory numbering and state management.
        
        The method performs the following steps:
        1. Loads simulation parameters from the original setup
        2. Creates Langevin integrator with original parameters
        3. Adds barostat for solvated systems (NPT ensemble)
        4. Loads the last saved state from XML file
        5. Sets up trajectory and log file reporters
        6. Runs simulation (time-based or step-based)
        7. Saves new state and updates restart configuration
        
        Side Effects:
            - Creates new trajectory file: output_traj_{last_state+1}.dcd
            - Appends to existing log file: log.txt
            - Saves updated state: last_state.xml
            - Saves final structure: last_state_{last_state+1}.pdb
            - Updates restart_setup.yml with incremented last_state counter
            - Terminates program execution with sys.exit(0)
            
        Raises:
            FileNotFoundError: If last_state.xml or restart_setup.yml files are missing.
            KeyError: If required parameters are missing from setup dictionary.
            
        Note:
            This method calls sys.exit(0) upon successful completion, terminating the program.
            The simulation continues exactly where the previous run left off, maintaining
            proper trajectory continuity and state consistency.
        """
        system = self.system
        modeller = self.modeller
        setup = self.setup 
       
        
        out_dir = self.config.restart
    
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

        if self.config.clock is not None:
            print('Starting simulation for', self.config.clock, ' mins ...')
            t1 = time.time()
            simulation.runForClockTime(self.config.clock * unit.minute)
            t2 = time.time()
            
        else :
            print('Starting simulation with', self.config.steps, 'steps ...')
            t1 = time.time()
            simulation.step(self.config.steps)
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



