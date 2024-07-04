#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem 8024
#SBATCH -o slurm.%N.%j.out      # STDOUT
#SBATCH -t 4:00:00              # time (D-HH:MM)
#SBATCH --gpus-per-node=1
#SBATCH --account=<def_somesposor>

module load StdEnv/2020
module load cuda/11.8
module load intel/2020.1.217

# number of restarts before exiting
max_restart=3
# restart number
restarts=0
# simulation time (for each restart)
sim_time=60
# simulation time  = sim_time * max_restart + sim_time 
# = 60 * 3 + 60 = 240 mins

source venv/bin/activate

python3 simulate.py -p test.pdb -l ATP --solvate --output ATP_explicit --clock $sim_time 

bash restart.sh ATP_explicit $sim_time $restarts $max_restart