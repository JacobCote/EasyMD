import sys, time, argparse
import os
import yaml
import pprint
import pickle
import openmm
from openmm.app import PDBFile, Simulation, StateDataReporter, DCDReporter
from openmm import app, unit, LangevinIntegrator 
import EasyMD.utils.utils as utils
#from prepare.prep_complex import prep_complex
#from prepare.prep_prot import prep_prot
#from restart.restart import prep_restart_ligand,restart_simulation,prep_restart
from EasyMD.utils.simulated_annealing import simulated_annealing
from EasyMD.simRunner import SimRunner
from EasyMD.sysGenerator import SysGenerator
from EasyMD.argManager import ArgManager
from .argManager.manager import ArgManager



def main():
    parser = argparse.ArgumentParser(
        prog='EasyMD',
        description="""EasyMD - Easy Molecular Dynamics Simulation Package
        
A comprehensive toolkit for setting up and running molecular dynamics simulations
with OpenMM. Supports protein-only and protein-ligand systems with both explicit
and implicit solvation.""",
        epilog="""Examples:
  # Basic protein simulation (10 ns, explicit solvent)
  python -m EasyMD --protein protein.pdb --steps 5000000 --solvate
  
  # Protein-ligand complex with custom parameters
  python -m EasyMD --protein complex.pdb --ligand ATP --steps 10000000 \\
                   --solvate --temperature 310 --ionic-strength 0.1
  
  # Fast implicit solvent simulation
  python -m EasyMD --protein protein.pdb --steps 1000000 --GBIS \\
                   --temperature 300 --step-size 0.002
  
  # Simulation with time-based duration
  python -m EasyMD --protein protein.pdb --clock 50ns --solvate
  
  # Keep crystal waters and use custom force fields
  python -m EasyMD --protein protein.pdb --steps 5000000 --solvate \\
                   --keep-water --protein-force-field amber99sb-ildn.xml
  
  # Restart from previous simulation
  python -m EasyMD --restart out_0/
  
  # Using configuration file
  python -m EasyMD --config my_simulation.yml

Output Files:
  - trajectory files: output_traj_*.dcd
  - log file: log.txt (energies, temperature, etc.)
  - final structure: last_state_*.pdb
  - system topology: topology.pkl
  - restart files: last_state.xml, restart_setup.yml

For more information and tutorials, visit: https://github.com/your-repo/EasyMD""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    argManager = ArgManager(parser)
    config = argManager.get_args()
    
    
    config.outdir = utils.get_outdir(config)
    

           



    # get the chosen or fastest platform
    platform = utils.get_platform()

    ## system preparation
    sysGenerator = SysGenerator(config)
    

    ## simulaton prep 
    sim = SimRunner(config,sysGenerator.modeller,sysGenerator.system)

    ## start sim
    sim.run()


if __name__ == "__main__":

    main()

