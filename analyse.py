import mdtraj as md
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import pickle

argparser = argparse.ArgumentParser(description='Analyse the trajectory')
argparser.add_argument('traj_directory', type=str, help='trajectory file')
argparser.add_argument('--rmsd', action='store_true', default=False, required=False, help='Calculate RMSD')
argparser.add_argument('--rmsf', action='store_true', default=False, required=False, help='Calculate RMSF')
args = argparser.parse_args()


## get number of dcd files in the directory

dir = args.traj_directory

files = [i for i in os.listdir(dir) if i.endswith('.dcd')]
n_dcd = len(files)

files = [ f'output_traj_{i}.dcd' for i in range(n_dcd) ]
print('Number of dcd files:', n_dcd)
print('Files:', files)
print('Directory:', dir)
print('Loading trajectories...')
## load topologies


print('Loading topologies...')
with open(dir + '/topology.pkl','rb') as f:
    topo_openmm = pickle.load(f)

print('Creating topology...')
topo = md.Topology.from_openmm(topo_openmm)
print('Done')
print('Loading trajectories....')
trajectories = [md.load(dir+'/'+i, top=topo) for i in files]
md_traj = md.join(trajectories)


if not os.path.isdir(dir+'/analysis'):
    os.mkdir(dir+'/analysis')



if args.rmsd:
    print('Calculating RMSD...')
    rmsd = md.rmsd(md_traj, md_traj, 0)
    plt.plot(rmsd)
    plt.title('RMSD')
    plt.xlabel('Frame')
    plt.ylabel('RMSD (Å^2)')
    plt.savefig(dir+'/analysis/rmsd.png')
    plt.clf()
    
if args.rmsf:
    print('Calculating RMSF...')
    
    # Calculate RMSF for the trajectory only CA atoms
    # select only CA atoms
    atoms_to_keep = [a.index for a in md_traj.topology.atoms if a.name == 'CA']
    
    rmsf = md.rmsf(md_traj ,md_traj, 0, atom_indices=atoms_to_keep )
    plt.plot(rmsf)
    plt.title('RMSF')
    plt.xlabel('Residue')
    plt.ylabel('RMSF (Å^2)')
    plt.savefig(dir+'/analysis/rmsf.png')
    plt.clf()
    
    
# get matrix of distances between centers of mass of all chains
if args.distances :
    import mdtraj as md
    
    pass