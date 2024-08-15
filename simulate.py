import sys, time, argparse
import os
import yaml
import pprint
import pickle
import openmm
from openmm.app import PDBFile, Simulation, StateDataReporter, DCDReporter
from openmm import app, unit, LangevinIntegrator 
import utils.utils as utils
from prepare.prep_complex import prep_complex
from prepare.prep_prot import prep_prot
from restart.restart import prep_restart_ligand,restart_simulation,prep_restart
from utils.simulated_annealing import simulated_annealing




t0 = time.time()


parser = argparse.ArgumentParser(description="Simulate", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-p", "--protein", required=False, help="Protein PDB file")
parser.add_argument("-l", "--ligand", required=False, help="Ligand name in pdb file (often LIG, check your pdb file to be sure of the name)")
parser.add_argument("-o", "--output", default=None, help="Name of an output directory")
parser.add_argument("-s", "--steps", type=int, default=None, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.001, help="Step size (ps")
parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("--solvate", action='store_true', help="Add solvent box")
parser.add_argument("--GBIS", action='store_true', help="Doesn't add solvent box, use Born generalize implicit solvent")
parser.add_argument("--padding", type=float, default=10, help="Padding for solvent box (A)")
parser.add_argument("--water-model", default="tip3p",
                    choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"],
                    help="Water model for solvation")
parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
parser.add_argument("--ionic-strength", type=float, default=0.1, help="Ionic strength for solvation")
parser.add_argument("--no-neutralize", action='store_true', help="Don't add ions to neutralize")
parser.add_argument("-e", "--equilibration-steps", type=int, default=200, help="Number of equilibration steps")
parser.add_argument("--protein-force-field", default='amber14-all.xml', help="Protein force field")
parser.add_argument("--ligand-force-field", default='openff-2.2.0', help="Ligand force field")
parser.add_argument("--water-force-field", default='amber/tip3p_standard.xml', help="Water force field")
parser.add_argument('--remove', nargs='*', help='Space separated molecules name to remove ex: --remove DMS LIG CA MG ... ', required=False, default=['DMS'])
parser.add_argument('--ph', type=float, help='Ph for the protonation state of the residus', required=False, default=7.0)
parser.add_argument("--restart", action='store_true', help="Use the program in restart mode.",default=False)
parser.add_argument("--restart_dir",type=str, help="path to the restart files", required=False, default='None')
parser.add_argument('--clock', type=float, help='Run the simulation based on clock time in minutes instead of steps.', required=False, default=None)
parser.add_argument('--simulated-annealing', action='store_true', help='Run a simulated annealing simulation', required=False, default=False)

args = parser.parse_args()
print("Simulate with these parameters: ")
pprint.pprint(vars(args))


# parser sanity check


if not args.simulated_annealing :

    if args.clock is not None and args.steps is not None:
        print('Please choose either --steps or --clock')
        exit(1)
    if args.steps is None and args.clock is None:
        print('Please provide either --steps or --clock')
        exit(1)

    if (args.solvate + args.GBIS % 2 == 0) and not args.restart :
        print('Please choose either --solvate or --GBIS')
        exit(1)
        
    if args.restart and  args.restart_dir is None:
        print('Please  directory with the restart files with --restart_dir <RESTARTDIR>')
        exit(1)
        
# get the chosen or fastest platform
platform = utils.get_platform()

# more specific forcefield kwargs, for now only the constraints
forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
last_state = 0

# restart from a previous simulation
if args.restart:    
    out_dir = args.restart_dir
    # load the setup file
    with open(out_dir+'/'+'restart_setup.yml') as f:
        setup = yaml.load(f, Loader=yaml.FullLoader)

    print('Restarting from state files in directory', args.restart_dir, 'with setup file : ' )
    pprint.pprint(setup)
    
    #check if ligand.sdf is present in the restart directory 
    if not os.path.isfile(out_dir+'/'+'ligand.sdf'):
        print('No ligand.sdf file found in the restart directory, preparing system without ligand ...')
        modeller, system = prep_restart(setup=setup, outdir=out_dir, forcefield_kwargs=forcefield_kwargs)
        
    else:
        print('Ligand.sdf file found in the restart directory, preparing system with ligand ...')
        modeller, system = prep_restart_ligand(setup=setup, outdir=out_dir, forcefield_kwargs=forcefield_kwargs)
    
    restart_simulation(
        system=system,
        modeller=modeller,
        restart_dir=args.restart_dir,
        setup=setup,
        clock=args.clock,
        step=args.steps,
         
    )
    
# create a new directory for the output if not restarting
else:
    if args.output is None:
   
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
        out_dir = args.output.rstrip('/')
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
           



pdb_in = args.protein
mol_in = args.ligand
output_traj_dcd = f'output_traj_{last_state}.dcd'
num_steps = args.steps
reporting_interval = args.interval
temperature = args.temperature * unit.kelvin
equilibration_steps = args.equilibration_steps

if args.ligand is None:
    modeller, system = prep_prot(
        pdb_in = pdb_in,
        list_of_molecules_to_remove=args.remove,
        solvate=args.solvate,
        protein_force_field=args.protein_force_field,
        water_force_field=args.water_force_field,
        water_model=args.water_model,
        positive_ion=args.positive_ion,
        negative_ion=args.negative_ion,
        ionic_strength=args.ionic_strength,
        no_neutralize=args.no_neutralize,
        padding=args.padding,
        ph=args.ph,
        outdir=out_dir,
        forcefield_kwargs=forcefield_kwargs
    )
    
else:
    modeller,system = prep_complex(
         pdb_in = pdb_in,
        list_of_molecules_to_remove=args.remove,
        lig_name=args.ligand,
        solvate=args.solvate,
        protein_force_field=args.protein_force_field,
        water_force_field=args.water_force_field,
        ligand_force_field=args.ligand_force_field,
        water_model=args.water_model,
        positive_ion=args.positive_ion,
        negative_ion=args.negative_ion,
        ionic_strength=args.ionic_strength,
        no_neutralize=args.no_neutralize,
        padding=args.padding,
        ph=args.ph,
        outdir=out_dir,
        forcefield_kwargs=forcefield_kwargs
    )


if args.simulated_annealing:
    simulated_annealing(modeller, system, temperature, out_dir, args.step_size, args.friction_coeff, args.interval, args.equilibration_steps)
    exit(0)                 
                            
# Set up the simulation
friction_coeff = args.friction_coeff / unit.picosecond
if args.clock is None:
    step_size = args.step_size * unit.picoseconds
    duration = (step_size * num_steps).value_in_unit(unit.nanoseconds)
    print('Simulating for {} ns'.format(duration))

integrator = LangevinIntegrator(temperature, friction_coeff, args.step_size * unit.picoseconds)

if args.solvate:
    system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature, 25))

if system.usesPeriodicBoundaryConditions():
    print('Default Periodic box: {}'.format(system.getDefaultPeriodicBoxVectors()))
else:
    print('No Periodic Box')
    


simulation = Simulation(modeller.topology, system, integrator)


context = simulation.context
context.setPositions(modeller.positions)



print('Minimising ...')
simulation.minimizeEnergy(maxIterations=10000)

# Write out the minimised PDB.
with open(out_dir+'/'+'minimised.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, context.getState(getPositions=True, ).getPositions(), file=outfile, keepIds=False)

# equilibrate
simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating ...')
simulation.step(equilibration_steps)

# Run the simulation.
# check for name 
if args.solvate:
    simulation.reporters.append(DCDReporter(out_dir+'/'+output_traj_dcd, reporting_interval, enforcePeriodicBox=False))
else :
    simulation.reporters.append(DCDReporter(out_dir+'/'+output_traj_dcd, reporting_interval))
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
#add a reporter for a log file
simulation.reporters.append(StateDataReporter(out_dir+'/'+'log.txt', reporting_interval, step=True, potentialEnergy=True, temperature=True))



# start simulation, either for a number of steps or for a clock time
if args.clock is not None:
    print('Starting simulation for', args.clock, 'mins ...')
    t1 = time.time()
    simulation.runForClockTime(args.clock * unit.minute)
    t2 = time.time()
else:
    print('Starting simulation with', num_steps, 'steps ...')
    t1 = time.time()
    simulation.step(num_steps)
    t2 = time.time()
    
# get duration of the simulation
step_s = args.step_size * unit.picoseconds
n_step = simulation.context.getStepCount()
duration = (n_step * step_s).value_in_unit(unit.nanoseconds)

#simulation.saveCheckpoint(out_dir+'/'+f'checkpoint_{last_state}.chk')
simulation.saveState(out_dir+'/'+'last_state.xml')
positions = simulation.context.getState(getPositions=True).getPositions()
# save the last state of the simulation
PDBFile.writeFile(simulation.topology, positions, open(out_dir+'/'+f'last_state_{last_state}.pdb', 'w'))
# save the last state 
PDBFile.writeFile(simulation.topology, positions, open(out_dir+'/'+'last_state.pdb', 'w'))
print('Simulation complete in {} mins at {}. Total wall clock time was {} mins'.format(
    round((t2 - t1) / 60, 3), temperature, round((t2 - t0) / 60, 3)))
print('Simulation time was', round(duration, 3), 'ns')

# save openmm topology
print('Writting topology file !!!')
filename = out_dir+'/topology.pkl'
with open(filename, 'wb')  as f :
    pickle.dump(simulation.topology,file = f)




print('Saving trajectory parameters in restart_setup.yml for eventual restart ...')
# write setup file for the next step
setup_dict = {
    'pdb': 'restart_model.pdb',
    'reporting_interval': args.interval,
    'step_size': args.step_size,
    'friction_coeff': args.friction_coeff,
    'temperature': args.temperature,
    'solvate': args.solvate,
    'state': out_dir+'/'+'last_state.xml',
    'protein_force_field': args.protein_force_field,
    'ligand_force_field': args.ligand_force_field,
    'water_force_field': args.water_force_field,
    'last_state': 0,
}



yaml.dump(setup_dict, open(out_dir+'/'+'restart_setup.yml', 'w'), default_flow_style=False)

#### POUR DEMAIN ####
# regarder pour le restart avec le state file et le restart_dir avec ou sans topo ?
