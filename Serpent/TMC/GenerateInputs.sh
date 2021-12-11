#!/bin/bash

#iso='U238'
#iso='H1'
#iso='Fe56'
#iso='O16'
#iso='U235'
#iso='Zr90'
iso='Ag107'

if [ "$iso" = "Ag107" ]
then
  nrand=99
else
  nrand=299
fi
nseed=4

#---
#  Generate SERPENT2 inputs (300 TMCs x 4 seeds x 3 control rods = 3600)
#---

# A single TMC sample is computed several times, with different seeds in order
# to have more statistics
for iseed in $(seq 1 4)
do
    for jiso in $(seq -w 000 $nrand)
	do
		isoid=$iso"_"$jiso
		for original in ../FullCore/TihangeFullCore*ppm.sss2
        do
            # Keep only the basename of the original file
            input=${original##*/}
            base=${input%.sss2}
			# Define a name for the new input file
            input=${input%.sss2}"_"$iseed"_"$isoid".sss2"
			cp $original $input
			# Replace an isotope with another TMC sample
            if [ "$iso" = "H1" ]
            then
			    sed -i "s/H1lwtr/$isoid/" $input
            else
			    sed -i "s/$iso/$isoid/" $input
            fi
            # Reduce to 100 inactive batches but use an already converged fission
            # source from a previous computation (unsampled, TMC-wise)
            sed -i "s/set pop.*/set pop 200000 2000 100 1.0/" $input
			echo 'src csw sf "../csw/'$base'.csw" -1' >> $input
            # By default, Serpent uses epoch as the RNG seed. For many calculations
            # launched at the same time, this implies many identical random seeds. We
            # add an integer to epoch, so that the various seeds are different.
            seed="set seed "$(date '+%s')$iseed
            echo $seed >> $input
        done
    done
done

#---
#  Generate corresponding SLURM instructions
#  They can be launched on Cobalt with
#  ls *_s*.sh | xargs -I '{}' ccc_msub '{}'
#---

for kslurm in {1..75}
do
    instructions=$iso'_s'$kslurm'.sh'
    # In each SLURM instruction, fully perform computations for 4 different TMC samples of that isotope
    jisomin=$(printf "%03d" $((4*($kslurm-1))))
    jisomax=$(printf "%03d" $((4*($kslurm)-1)))
    echo '#!/bin/sh                                                                        '  > $instructions
    # In each SLURM instruction, book two Skylake nodes (40 CPUs each)
    echo '#MSUB -n 8                                                                       ' >> $instructions
    echo '#MSUB -c 10                                                                      ' >> $instructions
    echo '#MSUB -q skylake                                                                 ' >> $instructions
    # On such a node, this computation should not take more than 20 hours
    echo '#MSUB -T 72000                                                                   ' >> $instructions
    echo '#MSUB -m scratch                                                                 ' >> $instructions
    echo '#MSUB -Q normal                                                                  ' >> $instructions
    echo 'cd ${BRIDGE_MSUB_PWD}                                                            ' >> $instructions
    echo 'sss2="../../../sss2_icc"                                                         ' >> $instructions
    echo 'iso="'$iso'"                                                                     ' >> $instructions
    echo 'for CB in 1206 1084 960                                                          ' >> $instructions
    echo 'do                                                                               ' >> $instructions
    echo '    for iseed in $(seq 1 4)                                                      ' >> $instructions
    echo '    do                                                                           ' >> $instructions
    echo '        for jiso in $(seq -w '$jisomin' '$jisomax')                              ' >> $instructions
    echo '        do                                                                       ' >> $instructions
    echo '            input="TihangeFullCore"$CB"ppm_"$iseed"_"$iso"_"$jiso".sss2"         ' >> $instructions
    echo '            ccc_mprun -E '"'"'--exclusive'"'"' -n 1 -c 10 $sss2 -omp 10 $input & ' >> $instructions
    echo '        done                                                                     ' >> $instructions
    echo '    done                                                                         ' >> $instructions
    echo '    wait  # Wait for all the ccc_mprun(s) to complete                            ' >> $instructions
    echo 'done                                                                             ' >> $instructions
    chmod 755 $instructions
done
