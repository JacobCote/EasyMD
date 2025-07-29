import sys
from openmm.app import PDBFile, Simulation, StateDataReporter, DCDReporter
from openmm import unit, LangevinIntegrator 


class AnnealingRunner:
    """
    Molecular dynamics simulation runner for simulated annealing protocols.
    
    This class implements simulated annealing, a technique that gradually increases
    temperature during simulation to help the system escape local energy minima
    and explore conformational space more effectively. This is particularly useful
    for protein folding studies, conformational sampling, and optimization problems.
    
    Attributes:
        config: Configuration object containing simulation parameters.
        modeller: OpenMM Modeller object with system topology and positions.
        system: OpenMM System object with forces and parameters.
        temperature (unit.Quantity): Starting simulation temperature in Kelvin.
        equilibration_steps (int): Number of equilibration steps.
        step_size (unit.Quantity): Integration time step in picoseconds.
    
    Example:
        >>> runner = AnnealingRunner(config, modeller, system)
        >>> runner.run()  # Starts the simulated annealing protocol
        
    Note:
        The annealing protocol gradually increases temperature from the base
        temperature by 0.1K increments over 2000 cycles of 100 steps each.
    """
    def __init__(self, config, modeller, system):
        """
        Initialize the AnnealingRunner with simulation parameters.
        
        Args:
            config: Configuration object containing simulation parameters including
                   temperature (starting point), intervals, and output settings.
            modeller: OpenMM Modeller object with system topology and initial positions.
            system: OpenMM System object containing force field parameters and constraints.
        """
        self.config = config
        self.modeller = modeller
        self.system = system
        self.temperature = self.config.temperature * unit.kelvin
        self.equilibration_steps = self.config.equilibration_steps
        self.step_size = self.config.step_size * unit.picoseconds
      

    def run(self):
        """
        Execute the simulated annealing molecular dynamics protocol.
        
        This method performs the following simulated annealing workflow:
        1. Sets up the Langevin integrator with initial temperature
        2. Performs energy minimization to remove bad contacts
        3. Equilibrates the system at starting temperature
        4. Runs annealing protocol: 2000 cycles of 100 steps each
        5. Gradually increases temperature by 0.1K per cycle
        6. Saves final annealed structure
        
        The annealing protocol helps the system escape local minima by providing
        thermal energy to overcome energy barriers, potentially finding lower
        energy conformations or exploring conformational space more thoroughly.
        
        Side Effects:
            - Creates trajectory file: simulated_anhealing.dcd
            - Creates log file: log.txt with energy and temperature data
            - Saves minimized structure: minimised.pdb
            - Saves final annealed structure: last_state.pdb
            - Terminates program execution with sys.exit(0)
            
        Raises:
            Exception: If energy minimization fails or simulation encounters errors.
            
        Note:
            The method calls sys.exit(0) upon completion, terminating the program.
            Temperature increases from base temperature to base + 200K over the protocol.
        """
        # Set up the simulation
        friction_coeff = self.config.friction_coeff / unit.picosecond
        
        #step_size = step_size * unit.picoseconds
        
    
        integrator = LangevinIntegrator(self.temperature, friction_coeff, self.step_size * unit.picoseconds)


        simulation = Simulation(self.modeller.topology, self.system, integrator)


        context = simulation.context
        context.setPositions(self.modeller.positions)



        print('Minimising ...')
        try:
            simulation.minimizeEnergy(maxIterations=10000)
            print('✓ Energy minimization completed successfully')
        except Exception as e:
            if "NaN" in str(e) or "coordinate is NaN" in str(e):
                print(f"❌ Energy minimization failed due to NaN coordinates: {e}")
                print("This usually indicates:")
                print("  - Overlapping atoms in the initial structure")
                print("  - Invalid ligand coordinates")
                print("  - Issues with the force field parameterization")
                print("\nTroubleshooting suggestions:")
                print("  1. Check your input PDB structure for overlapping atoms")
                print("  2. Verify ligand coordinates are reasonable")
                print("  3. Try using --keep-water flag if crystal waters are important")
                print("  4. Consider using a different force field")
                raise RuntimeError(f"Energy minimization failed due to coordinate issues: {e}")
            else:
                print(f"❌ Energy minimization failed: {e}")
                raise

        # Write out the minimised PDB.
        with open(self.config.out_dir+'/'+'minimised.pdb', 'w') as outfile:
            PDBFile.writeFile(self.modeller.topology, context.getState(getPositions = True, ).getPositions(), file=outfile, keepIds=False)

        # equilibrate
        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Equilibrating ...')
        simulation.step(self.equilibration_steps)

        # Run the simulation.
        # check for name 
    
        
        simulation.reporters.append(DCDReporter(self.out_dir+'/'+'simulated_anhealing.dcd', self.config.interval))
        simulation.reporters.append(StateDataReporter(sys.stdout, self.config.interval * 5, step=True, potentialEnergy=True, temperature=True))
        #add a reporter for a log file
        simulation.reporters.append(StateDataReporter(self.config.out_dir+'/'+'log.txt', self.config.interval, step=True, potentialEnergy=True, temperature=True))


        ## simulate 
        for i in range(2000):
            
            integrator.setTemperature(self.temperature+(0.1*i)* unit.kelvin)
            simulation.step(100)


        # get duration of the simulation
        step_s = self.step_size * unit.picoseconds
        n_step = simulation.context.getStepCount()
        duration = (n_step * step_s).value_in_unit(unit.nanoseconds)
        positions = simulation.context.getState(getPositions=True).getPositions()

        #simulation.saveCheckpoint(out_dir+'/'+f'checkpoint_{last_state}.chk')
        
        # save the last state of the simulation
        PDBFile.writeFile(simulation.topology, positions, open(self.config.out_dir+'/'+f'last_state.pdb', 'w'))
        exit(0)
        
        
    