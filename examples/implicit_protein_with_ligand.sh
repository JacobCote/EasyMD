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

# you may need to conda init if on a cluster
conda init
conda activate mdEnv
# number of restarts before exiting
max_restart=3
# restart number
restart=0
# simulation time (for each restart)
sim_time=200
# simulation time  = sim_time * max_restart + sim_time. This setup is made for running on a cluster with a time limit of 4 hours
#  200 * 3 + 200 = 800 mins of calculation time split into 4 ish hours block

source venv/bin/activate

python3 simulate.py -p test.pdb -l ATP --GBIS --output ATP_implicit --clock $sim_time --remove MG

# restart the simulation if restarts is less than max_restart
if [ $restarts -lt $max_restart ]; then
    bash restart.sh ATP_explicit $sim_time $(($restarts+1)) $max_restart
fi
