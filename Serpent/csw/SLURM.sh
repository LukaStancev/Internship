#!/bin/bash

sss2=/soft_snc/SERPENT2/2.1.32/sss2
if ! command -v $sss2 &> /dev/null
then
    echo "Serpent 2 could not be found in $sss2"
    exit 1
fi

for original in ../FullCore/TihangeFullCore*ppm.sss2
do
    # Keep only the basename of the original file
    input=${original##*/}
    cp $original $input
    csw=${input%.sss2}".csw"
    echo "set csw "$csw >> $input
    sbatch -c 40 --partition="seq,dev,par_IB" --mem=40G --wrap="srun $sss2 -omp 40 $input"
done
