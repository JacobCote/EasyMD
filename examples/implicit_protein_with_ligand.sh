#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem 4012
#SBATCH -o slurm.%N.%j.out      # STDOUT
#SBATCH -t 4:00:00              # time (D-HH:MM)
#SBATCH --gpus-per-node=1
#SBATCH --account=<def_somesposor>


# you need to move this script and restart.sh to the root directory of the project

module load StdEnv/2020
module load cuda/11.8
module load intel/2020.1.217

# you may need to conda init if on a cluster
conda init bash
conda activate mdEnv
# number of restarts before exiting
max_restart=3
# restart number
restarts=0
# simulation time (for each restart)
sim_time=200

outdir=test_implicit
# simulation time  = sim_time * max_restart + sim_time. This setup is made for running on a cluster with a time limit of 4 hours
#  200 * 3 + 200 = 800 mins of calculation time split into 4 ish hours block

# here we also remove MG ions from the simulation because of implicit solvent
python3 simulate.py -p test.pdb -l GNP --GBIS --output $outdir --clock $sim_time --remove MG

# restart the simulation if restarts is less than max_restart
if [ $restarts -lt $max_restart ]; then
    bash restart.sh $outdir $sim_time $(($restarts+1)) $max_restart
fi
