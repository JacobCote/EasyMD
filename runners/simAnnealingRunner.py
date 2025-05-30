import sys
from openmm.app import PDBFile, Simulation, StateDataReporter, DCDReporter
from openmm import unit, LangevinIntegrator 


class AnnealingRunner():
    def __init__(self,config,modeller,system):
        self.config = config
        self.modeller = modeller
        self.system = system
        self.temperature = self.config.temperature * unit.kelvin
        self.equilibration_steps = self.config.equilibration_steps
        self.step_size = self.config.step_size * unit.picoseconds
      


    def simulated_annealing(self,):
                    
        # Set up the simulation
        friction_coeff = friction_coeff / unit.picosecond
        
        #step_size = step_size * unit.picoseconds
        
    
        integrator = LangevinIntegrator(self.temperature, friction_coeff, self.step_size * unit.picoseconds)


        simulation = Simulation(self.modeller.topology, self.system, integrator)


        context = simulation.context
        context.setPositions(self.modeller.positions)



        print('Minimising ...')
        simulation.minimizeEnergy(maxIterations=10000)

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
        
        
    