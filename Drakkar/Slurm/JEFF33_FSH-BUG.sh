#!/bin/bash
#SBATCH -c 1  ## Number of cores to be reserved
#SBATCH -t 0-06:00:00
#SBATCH --mem-per-cpu=4G

#rm -rf ../Linux_x86_64
cd ../data

srun ../runVersion5.sh Fessenheim-Bugey.x2m

