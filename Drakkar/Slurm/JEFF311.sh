#!/bin/bash
#SBATCH -c 1  ## Number of cores to be reserved
#SBATCH -t 0-06:00:00
#SBATCH --mem-per-cpu=4G

rm -rf Linux_x86_64
cd data

srun ../runVersion5.sh JEFF311.x2m

grep --color=none "RELATIVE POWER" -A 16 ../Linux_x86_64/JEFF311.result | grep --color=none "  1 " -A 14
grep --color=none "EFFECTIVE M" ../Linux_x86_64/JEFF311.result | cut -d"=" -f2
