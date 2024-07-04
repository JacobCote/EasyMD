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


def prep_restart_ligand(setup: dict,outdir : str, forcefield_kwargs: dict ) -> Tuple[Modeller, SystemGenerator]:
    '''
    prepare system for restart simulation with ligand
    :param setup: dict, setup parameters for the simulation wich is created when the simulation is started for the first time
    :param outdir: str, output directory where the restart files are stored
    :param forcefield_kwargs: dict, forcefield parameters
    :return: Tuple[Modeller, SystemGenerator], Modeller object and SystemGenerator object
    '''
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




def prep_restart(setup: dict,outdir : str, forcefield_kwargs: dict ) -> Tuple[Modeller, SystemGenerator]:
    '''
    prepare system for restart simulation withoput ligand
    :param setup: dict, setup parameters for the simulation wich is created when the simulation is started for the first time
    :param outdir: str, output directory where the restart files are stored
    :param forcefield_kwargs: dict, forcefield parameters
    :return: Tuple[Modeller, SystemGenerator], Modeller object and SystemGenerator object
    '''
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



def restart_simulation(system, modeller,restart_dir, setup,  clock=None, step=None):
    '''
    restart the simulation from the last state of the previous simulation from a stat.xml file and topology
    :param system: openmm.System, the system object
    :param modeller: openmm.Modeller, the modeller object
    :param restart_dir: str, the directory where the restart files are stored
    :param setup: dict, setup parameters for the simulation wich is created when the simulation is started for the first time
    :param clock: int, the time in minutes for the simulation to run
    :param step: int, the number of steps for the simulation to run
    :return: None
    '''
    
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



