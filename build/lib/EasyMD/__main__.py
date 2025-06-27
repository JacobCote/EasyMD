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


    #t0 = time.time()

    # setup argManager
    parser = argparse.ArgumentParser(description="Simulate", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argManager = ArgManager(parser)
    config = argManager.getargs()
    print("Simulate with these parameters: ")
    pprint.pprint(vars(config))

    # get the chosen or fastest platform
    platform = utils.get_platform()


    ## system preparation
    modeller, system = SysGenerator(config)

    ## simulaton prep 
    sim = SimRunner(config,modeller,system)

    ## start sim
    sim.run()


if __name__ == "__main__":
    print("ok")

    main()

