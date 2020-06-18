#!/bin/bash
#SBATCH -c 1  ## Number of cores to be reserved
#SBATCH -t 0-06:00:00
#SBATCH --mem-per-cpu=4G

rm -rf Reflector/Linux_x86_64 Reflector/save/*
cd Reflector/data

srun dragon TousPaliers.x2m v5bev1761
