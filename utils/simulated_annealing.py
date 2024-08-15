import sys
from openmm.app import PDBFile, Simulation, StateDataReporter, DCDReporter
from openmm import app, unit, LangevinIntegrator 





def simulated_annealing(modeller, system, temperature, out_dir,step_size, friction_coeff, reporting_interval, equilibration_steps):
                
    # Set up the simulation
    friction_coeff = friction_coeff / unit.picosecond
    
    #step_size = step_size * unit.picoseconds
    
   
    integrator = LangevinIntegrator(temperature, friction_coeff, step_size * unit.picoseconds)


    simulation = Simulation(modeller.topology, system, integrator)


    context = simulation.context
    context.setPositions(modeller.positions)



    print('Minimising ...')
    simulation.minimizeEnergy(maxIterations=10000)

    # Write out the minimised PDB.
    with open(out_dir+'/'+'minimised.pdb', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, context.getState(getPositions = True, ).getPositions(), file=outfile, keepIds=False)

    # equilibrate
    simulation.context.setVelocitiesToTemperature(temperature)
    print('Equilibrating ...')
    simulation.step(equilibration_steps)

    # Run the simulation.
    # check for name 
   
    
    simulation.reporters.append(DCDReporter(out_dir+'/'+'simulated_anhealing.dcd', reporting_interval))
    simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
    #add a reporter for a log file
    simulation.reporters.append(StateDataReporter(out_dir+'/'+'log.txt', reporting_interval, step=True, potentialEnergy=True, temperature=True))


    ## simulate 
    for i in range(2000):
        
        integrator.setTemperature(temperature+(0.1*i)* unit.kelvin)
        simulation.step(100)


    # get duration of the simulation
    step_s = step_size * unit.picoseconds
    n_step = simulation.context.getStepCount()
    duration = (n_step * step_s).value_in_unit(unit.nanoseconds)
    positions = simulation.context.getState(getPositions=True).getPositions()

    #simulation.saveCheckpoint(out_dir+'/'+f'checkpoint_{last_state}.chk')
    
    # save the last state of the simulation
    PDBFile.writeFile(simulation.topology, positions, open(out_dir+'/'+f'last_state.pdb', 'w'))
    exit(0)
    
    
  
