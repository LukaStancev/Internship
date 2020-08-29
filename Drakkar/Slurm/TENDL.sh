#!/bin/bash
# Usage : ./TENDL.sh
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

intg=0
for iso in O16 U235 U238 Zr90 Zr91 Zr94 Zr96
do
  if [ "$iso" = "Zr91" ]
  then
    maxrand=180
  elif [ "$iso" = "Zr92" ] || [ "$iso" = "Zr94" ]
  then
    maxrand=200
  elif [ "$iso" = "Zr92" ]
  then
    maxrand=250
  else
    maxrand=299
  fi
  for (( irand=0; irand<= $maxrand; irand++ ))
  do
    # Unique integer to avoid collisions in runVersion5.sh
    intg=$((intg+1))
    # Prepare the files to be modified
    cp ../data/Tihange.x2m    Tihange_${iso}_${irand}.x2m
    cp ../data/Tihange.access Tihange_${iso}_${irand}.access
    cp ../data/Tihange.save   Tihange_${iso}_${irand}.save
    # Change the isotope from JEFF-3.1.1 to a random one from TENDL
    sed -i 's/'$iso'.*/'$iso' '$irand' ;/' Tihange_${iso}_${irand}.x2m
    # Writing Slurm instruction file
    echo "#!/bin/bash"                                                 > ${iso}_${irand}_Tihange.sh
    echo "#SBATCH -c1  ## Number of cores to be reserved"             >> ${iso}_${irand}_Tihange.sh
    echo "#SBATCH -t 0-00:30:00"                                      >> ${iso}_${irand}_Tihange.sh
    echo "#SBATCH --mem-per-cpu=4G"                                   >> ${iso}_${irand}_Tihange.sh
    echo "srun ../runVersion5.sh Tihange_${iso}_${irand}.x2m ${intg}" >> ${iso}_${irand}_Tihange.sh
    # Executing the instruction file
    chmod 755 ${iso}_${irand}_Tihange.sh
    sbatch ${iso}_${irand}_Tihange.sh
  done
done
