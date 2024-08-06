#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem 4012
#SBATCH -o slurm.%N.%j.out      # STDOUT
#SBATCH -t 4:00:00              # time (D-HH:MM)
#SBATCH --gpus-per-node=1
#SBATCH --account=<def_somesposor>



module load StdEnv/2020
module load cuda/11.8
module load intel/2020.1.217

# you may need to conda init if on a cluster
conda init bash
conda activate mdEnv

restart_dir=$1
clock=$2
restarts=$3
max_restart=$4


source venv/bin/activate
python3 simulate.py --restart --restart_dir $restart_dir --clock $clock

# if restarts is less than max_restart, restart the simulation 
if [ $restarts -lt $max_restart ]; then
    bash restart.sh $restart_dir $clock $(($restarts+1)) $max_restart
fi