#!/bin/sh

sss2=/soft_snc/SERPENT2/2.1.32/sss2
if ! command -v $sss2 &> /dev/null
then
    echo "Serpent 2 could not be found in $sss2"
    exit 1
fi

for input in TihangeFullCore*.sss2
do
    echo $input
    sbatch -c 8 --mem-per-cpu=8G --wrap="srun $sss2 -omp 8 $input"
done
