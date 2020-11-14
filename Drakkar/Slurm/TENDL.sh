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
for iso in O16 U235 U238 Zr90 Zr91 Zr92 Zr94 Zr96 Ni58 Fe54 In115 Ag107 Ag109 Cd106 Cd108 Cd110 Cd111 Cd112 Cd113 Cd114 Cd116 Fe56 Cr52
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
  elif [ "$iso" == 'Fe54' ]
  then
    maxrand=180
  elif [ "$iso" == 'In115' ] || [ "$iso" == 'Ag107' ]
  then
    maxrand=100
  elif [ "$iso" == 'Ag109' ] || [ "$iso" == 'Cd112' ] || [ "$iso" == 'Cd114' ]
  then
    maxrand=60
  elif [ "$iso" == 'Cd110' ] || [ "$iso" == 'Cd111' ] || [ "$iso" == 'Cd113' ]
  then
    maxrand=40
  elif [ "$iso" == 'Cd106' ] || [ "$iso" == 'Cd108' ] || [ "$iso" == 'Cd116' ]
  then
    maxrand=20
  elif [ "$iso" == 'Cr52' ] || [ "$iso" == 'Fe56' ]
  then
    maxrand=10
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
