#!/bin/bash
# Usage : FSH-BUG_TMC_with_interactions.sh

nrand=2000

# List of isotopes that are being randomly sampled
isotopes="H1_H2O B10 B11 O16 U235 U238 Zr90 Zr91 Zr92 Zr94 Zr96 Ni58 Fe54 Fe56 Cr52"

# Can be either 'Tihange' or 'Fessenheim-Bugey'
case=Fessenheim-Bugey

# Function generating pseudorandom integers from Python (i.e. Mersenne Twister).
#
# Usage : rand seed maxrange nrand
# where
# * seed     : the number used to initialize the pseudorandom number generator.
# * maxrange : the maximum integer that can be produced (inclusive). The minimum integer is zero (also included).
# * nrand    : the number of pseudorandom numbers to be produced.
# Ref. : https://docs.python.org/fr/3.10/library/random.html
rand(){
    echo $(python -c "import random;random.seed($1);print([random.randrange($2) for i in range(0,$3)])")
}

# Prepare an empty directory
rm -rf ../Linux_x86_64
inum=1
while [ -d ../randomized_data_$inum ]
  do
  inum=`expr $inum + 1 `
done
LaunchDir=../randomized_data_$inum
mkdir $LaunchDir
echo "Launch Directory:" $LaunchDir
cd $LaunchDir

cp ../data/isotconvlist.csv .

# nrand identical files are produced, initially bearing -33 for all isotopes (=JEFF-3.3 best estimate)
intg=0
for (( i=0; i<= $nrand-1; i++ )) # Until nrand-1 because 0 is included. Total is indeed nrand.
do
    # Prepare the files to be modified
    cp ../data/${case}.x2m    ${case}_${i}.x2m
    cp ../data/${case}.access ${case}_${i}.access
    cp ../data/${case}.save   ${case}_${i}.save
    # Unique integer to avoid collisions in runVersion5.sh
    intg=$((intg+1))
    # Writing Slurm instruction file, but not launching it
    echo "#!/bin/bash"                                                 > ${i}_${case}.sh
    echo "#SBATCH -c1  ## Number of cores to be reserved"             >> ${i}_${case}.sh
    echo "#SBATCH -t 0-00:30:00"                                      >> ${i}_${case}.sh
    echo "#SBATCH --mem-per-cpu=4G"                                   >> ${i}_${case}.sh
    echo "srun ../runVersion5.sh ${case}_${i}.x2m ${intg}" >> ${i}_${case}.sh
    chmod 755 ${i}_${case}.sh
done

# The use of a predefined seed allows reproducibility
seed=1354

for iso in $isotopes
do
  if [ "$iso" = "Zr91" ] || [ "$iso" == 'Fe54' ]
  then
    maxrand=179
  elif [ "$iso" = "Zr92" ] || [ "$iso" = "Zr94" ] || [ "$iso" = "Zr96" ]
  then
    maxrand=199
  elif [ "$iso" == 'Cr52' ]
  then
    maxrand=9
  else
    maxrand=299
  fi
  echo $iso
  # For each isotope, nrand pseudorandom numbers are generated (between 0 and maxrand) and injected into each of the nrand files
  rands=$(rand $seed $maxrand $nrand)
  i=0
  for rawrand in $rands
  do
      # Filter out anything given by Python that's not a number (commas, brackets...)
      rand=$(echo $rawrand | sed 's/[^0-9]*//g')
      # Change the isotope evaluation from JEFF-3.3 best-estimate to a randomly sampled one
      sed -i 's/'$iso'.*/'$iso' '$rand' ;/' ${case}_${i}.x2m
      # The random seed is updated in order to avoid producing the same numbers every time
      seed=$(($seed + $rand))
      i=$((i+1))
  done

done

# Launching computations
for (( i=0; i<= $nrand-1; i++ ))
do
    sbatch ${i}_${case}.sh
done
