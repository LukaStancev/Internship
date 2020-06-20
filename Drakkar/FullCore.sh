#!/bin/bash
#SBATCH -c 1  ## Number of cores to be reserved
#SBATCH -t 0-06:00:00
#SBATCH --mem-per-cpu=4G

rm -rf FullCore/Linux_x86_64
cd FullCore/data

srun donjon Tihange.x2m v5bev1879

grep "RELATIVE POWER" -A 16 ../Linux_x86_64/Tihange.result | grep "  1 " -A 14
grep "EFFECTIVE M" ../Linux_x86_64/Tihange.result | cut -d"=" -f2
