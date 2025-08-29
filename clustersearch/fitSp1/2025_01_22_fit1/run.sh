#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH --array=1-400:1
#SBATCH -t 00-03:00                         # Runtime in D-HH:MM format
#SBATCH --qos=short                           # Partition to run in
#SBATCH --mem=100                    # Memory total in MB (for all cores)
#SBATCH -o %A_%a.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e %A_%a.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=NONE   # Email to which notifications will be sent
#SBATCH -J SP1fit


module load Python/3.12.4
source /users/romartinez/romartinez/venvs/pythonfits/bin/activate 

#python makescratchdir.py
python fit.py $SLURM_ARRAY_TASK_ID
