#!/bin/bash

sss2=/soft_snc/SERPENT2/2.1.32/sss2
if ! command -v $sss2 &> /dev/null
then
    echo "Serpent 2 could not be found in $sss2"
    exit 1
fi

# 50 random samples = 10% uncertainty (1 sigma) of the standard deviation
for i in {1..50}
do
    for original in ../FullCore/TihangeFullCore*ppm.sss2
    do
        # Keep only the basename of the original file
        input=${original##*/}
        # Define a name for the new input file
        input=${input%.sss2}"_"$i".sss2"
        cp $original $input
        # Adjust the options in order to collect the detector results
        # obtained at each generation
        sed -i 's/set pop.*/set pop 100000 1200 1 1.0/' $input
        sed -i 's/^det .*/& dhis 1/' $input
        echo "set his 1" >> $input
        # By default, Serpent uses epoch as the RNG seed. For many
        # calculations launched at the same time, this implies many
        # identical random seeds. We add an integer to epoch, so that the
        # various seeds are different.
        seed="set seed "$(date '+%s')$i
        echo $seed >> $input
    done
done

for input in TihangeFullCore*ppm_*.sss2
do
    echo $input
    sbatch -c 5 --mem=30G --partition="seq,dev,par_IB" --exclude="../Exclude.txt" --wrap="srun $sss2 -omp 5 $input"
done
