import mdtraj as md
import argparse
import os
import numpy as np 
import pandas as pd
import pickle

argparser = argparse.ArgumentParser(description='Extract CA atoms from the trajectory')
argparser.add_argument('traj_directory', type=str, help='trajectory file')
args = argparser.parse_args()
# create a h5 dir
dir = args.traj_directory
if not os.path.exists(dir + '/CA-only'):
    os.mkdir(dir + '/CA-only')



files = [i for i in os.listdir(dir) if i.endswith('.dcd')]

files = [ f'output_traj_{i}.dcd' for i in range(len(files)) ]

print('Files:', files)

print('Loading topologies...')
with open(dir + '/topology.pkl','rb') as f:
    topo_openmm = pickle.load(f)

print('Creating topology...')
topo = md.Topology.from_openmm(topo_openmm)
print('Done')

for i, f in enumerate(files):
    print(f'processing file number {i} out of {len(files)}')
    print('loading trajectory...')
    traj = md.load(dir+'/'+f, top=topo)
    print('Extracting CA ...')
        
    res_2_keep = set([a.index for a in traj.topology.residues if not a.is_water ])

    atoms_to_keep = [a.index for a in traj.topology.atoms if a.residue.index in res_2_keep  ]


    traj.restrict_atoms(atoms_to_keep)
    # this acts inplace on the trajectory

    traj.save(dir + f'/CA-only/CA-only_{i}.dcd')
    traj[0].save(dir + '/CA-only/CA-Only.pdb')
    del traj

    