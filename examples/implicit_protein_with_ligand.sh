#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem 8024
#SBATCH -o slurm.%N.%j.out      # STDOUT
#SBATCH -t 4:00:00              # time (D-HH:MM)
#SBATCH --gpus-per-node=1
#SBATCH --account=<def_somesposor>


module load StdEnv/2020
module load cuda/11.4
module load intel/2020.1.217
# number of restarts before exiting
max_restart=3
# restart number
restart=0
# simulation time (for each restart)
sim_time=60
# simulation time  = sim_time * max_restart + sim_time 
# = 60 * 3 + 60 = 240

source venv/bin/activate

python3 simulate.py -p test.pdb -l ATP --GBIS --output ATP_implicit --clock $sim_time 

bash restart.sh ATP_imnplicit $sim_time $restart $max_restart