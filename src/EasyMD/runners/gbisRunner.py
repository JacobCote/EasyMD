import sys
from openmm.app import PDBFile, Simulation, StateDataReporter, DCDReporter
from openmm import unit, LangevinIntegrator 
import openmm
import pickle
import time
import yaml

class GBISRunner:
    """
    Molecular dynamics simulation runner for implicit solvent (GBIS) systems.
    
    This class handles MD simulations of protein or protein-ligand complexes using
    the Generalized Born Implicit Solvent (GBIS) model. It provides faster simulations
    compared to explicit solvent by treating solvent effects implicitly, making it
    suitable for conformational sampling and preliminary studies.
    
    Attributes:
        config: Configuration object containing simulation parameters.
        modeller: OpenMM Modeller object with system topology and positions.
        system: OpenMM System object with forces and parameters.
        temperature (unit.Quantity): Simulation temperature in Kelvin.
        equilibration_steps (int): Number of equilibration steps.
        step_size (unit.Quantity): Integration time step in picoseconds.
    
    Example:
        >>> runner = GBISRunner(config, modeller, system)
        >>> runner.run()  # Starts the complete GBIS simulation workflow
        
    Note:
        GBIS simulations do not use periodic boundary conditions and do not require
        a barostat since there is no explicit solvent box.
    """

    def __init__(self, config, modeller, system):
        """
        Initialize the GBISRunner with simulation parameters.
        
        Args:
            config: Configuration object containing simulation parameters including
                   temperature, steps, intervals, force fields, and output settings.
            modeller: OpenMM Modeller object with system topology and initial positions.
            system: OpenMM System object containing GBIS force field parameters and constraints.
        """
        self.config = config
        self.modeller = modeller
        self.system = system
        self.temperature = self.config.temperature * unit.kelvin
        self.equilibration_steps = self.config.equilibration_steps
        self.step_size = self.config.step_size * unit.picoseconds
      

    def run(self):
        """
        Execute the complete molecular dynamics simulation workflow for GBIS systems.
        
        This method performs the following steps for implicit solvent simulations:
        1. Sets up the Langevin integrator with temperature and friction control
        2. Performs energy minimization (fewer iterations than explicit solvent)
        3. Equilibrates the system at target temperature
        4. Runs production simulation (time-based or step-based)
        5. Saves final state, trajectory, and restart configuration
        
        Unlike explicit solvent simulations, GBIS does not require:
        - Monte Carlo barostat (no pressure control needed)
        - Periodic boundary conditions
        - Extensive energy minimization
        
        Side Effects:
            - Creates trajectory file: output_traj_0.dcd
            - Creates log file: log.txt with energy and temperature data
            - Saves minimized structure: minimised.pdb
            - Saves final state: last_state.xml and last_state_0.pdb
            - Saves topology: topology.pkl (pickled OpenMM topology)
            - Creates restart configuration: restart_setup.yml
            
        Raises:
            Exception: If energy minimization fails or simulation encounters errors.
            
        Note:
            GBIS simulations typically run faster than explicit solvent but may be
            less accurate for systems where specific solvent interactions are important.
        """
        # Set up the simulation
        last_state = 0
        output_traj_dcd = f'output_traj_{last_state}.dcd'
        friction_coeff = self.config.friction_coeff / unit.picosecond
        if self.config.clock is None:
            step_size = self.config.step_size * unit.picoseconds
            duration = (step_size * self.config.steps).value_in_unit(unit.nanoseconds)
            print('Simulating for {} ns'.format(duration))

        integrator = LangevinIntegrator(self.temperature, friction_coeff, self.config.step_size * unit.picoseconds)

        if self.system.usesPeriodicBoundaryConditions():
            print('Default Periodic box: {}'.format(self.system.getDefaultPeriodicBoxVectors()))
        else:
            print('No Periodic Box')
            


        simulation = Simulation(self.modeller.topology, self.system, integrator)


        context = simulation.context
        context.setPositions(self.modeller.positions)



        print('Minimising ...')
        simulation.minimizeEnergy(maxIterations=10)

        # Write out the minimised PDB.
        with open(self.config.outdir+'/'+'minimised.pdb', 'w') as outfile:
            PDBFile.writeFile(self.modeller.topology, context.getState(getPositions=True, ).getPositions(), file=outfile, keepIds=False)

        # equilibrate
        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Equilibrating ...')
        simulation.step(self.config.equilibration_steps)

        # Run the simulation.
   
        simulation.reporters.append(DCDReporter(self.config.outdir+'/'+output_traj_dcd, self.config.interval))
        simulation.reporters.append(StateDataReporter(sys.stdout, self.config.interval * 5, step=True, potentialEnergy=True, temperature=True))
        #add a reporter for a log file
        simulation.reporters.append(StateDataReporter(self.config.outdir+'/'+'log.txt', self.config.interval, step=True, potentialEnergy=True, temperature=True))



        # start simulation, either for a number of steps or for a clock time
        if self.config.clock is not None:
            print('Starting simulation for', self.config.clock, 'mins ...')
            t1 = time.time()
            simulation.runForClockTime(self.config.clock * unit.minute)
            t2 = time.time()
        else:
            print('Starting simulation with', self.config.steps, 'steps ...')
            t1 = time.time()
            simulation.step(self.config.steps)
            t2 = time.time()
            
        # get duration of the simulation
        step_s = self.config.step_size * unit.picoseconds
        n_step = simulation.context.getStepCount()
        duration = (n_step * step_s).value_in_unit(unit.nanoseconds)

        #simulation.saveCheckpoint(self.config.outdir+'/'+f'checkpoint_{last_state}.chk')
        simulation.saveState(self.config.outdir+'/'+'last_state.xml')
        positions = simulation.context.getState(getPositions=True).getPositions()
        # save the last state of the simulation
        PDBFile.writeFile(simulation.topology, positions, open(self.config.outdir+'/'+f'last_state_{last_state}.pdb', 'w'))
        # save the last state 
        PDBFile.writeFile(simulation.topology, positions, open(self.config.outdir+'/'+'last_state.pdb', 'w'))
        print('Simulation complete in {} mins at {}. '.format(
            round((t2 - t1) / 60, 3), self.temperature,))
        print('Simulation time was', round(duration, 3), 'ns')

        # save openmm topology
        print('Writting topology file !!!')
        filename = self.config.outdir+'/topology.pkl'
        with open(filename, 'wb')  as f :
            pickle.dump(simulation.topology,file = f)




        print('Saving trajectory parameters in restart_setup.yml for eventual restart ...')
        # write setup file for the next step
        setup_dict = {
            'pdb': 'restart_model.pdb',
            'reporting_interval': self.config.interval,
            'step_size': self.config.step_size,
            'friction_coeff': self.config.friction_coeff,
            'temperature': self.config.temperature,
            'solvate': self.config.solvate,
            'state': self.config.outdir+'/'+'last_state.xml',
            'protein_force_field': self.config.protein_force_field,
            'ligand_force_field': self.config.ligand_force_field,
            'water_force_field': self.config.water_force_field,
            'last_state': 0,
        }



        yaml.dump(setup_dict, open(self.config.outdir+'/'+'restart_setup.yml', 'w'), default_flow_style=False)

        
        
    