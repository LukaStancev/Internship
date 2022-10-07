#!/bin/sh

sss2=/soft_snc/SERPENT2/2.1.32/sss2
if ! command -v $sss2 &> /dev/null
then
    echo "Serpent 2 could not be found in $sss2"
    exit 1
fi

for c in {40..1}
do
    echo "c="$c
    for ppm in 1206 1084 960
    do
        echo "ppm="$ppm
        input="TihangeFullCore"$ppm"ppm_"$c".sss2"
        cp "TihangeFullCore"$ppm"ppm.sss2" $input
        srun -c $c --mem=40G --partition="seq,dev,par_IB" --exclude="./node119.txt" $sss2 -omp $c $input
    done
done
