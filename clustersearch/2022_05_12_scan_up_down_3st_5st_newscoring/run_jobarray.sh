#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH --array 1-221:1
#SBATCH -t 30-00:00                         # Runtime in D-HH:MM format
#SBATCH -p long                           # Partition to run in
#SBATCH --mem=120                          # Memory total in MB (for all cores)
#SBATCH -o %A_%a.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e %A_%a.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=abc123@hms.harvard.edu   # Email to which notifications will be sent
#SBATCH -J du_pc 

module load gcc/6.2.0
module load python/3.6.0
source /home/rm335/myenvpy36/bin/activate

#python makescratchdir.py
python boundaries_scan1.py $SLURM_ARRAY_TASK_ID
