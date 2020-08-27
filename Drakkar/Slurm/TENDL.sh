#!/bin/bash
# Usage : ./TENDL.sh
rm -rf ../Linux_x86_64
CurrentDir=$PWD
inum=1
while [ -d /tmp/launchdir$inum ]
  do
  inum=`expr $inum + 1 `
done
LaunchDir=/tmp/launchdir$inum
mkdir $LaunchDir
echo "Launch Directory:" $LaunchDir
cd $LaunchDir

ln -s $CurrentDir/../draglib* .
mkdir proc
cp $CurrentDir/../proc/* proc/.
mkdir data
cd data
cp $CurrentDir/../data/isotconvlist.csv .

for iso in O16 U235 U238 Zr90 Zr91 Zr94 Zr96
do
  if [ "$iso" -eq "Zr91" ]
  then
    maxrand=180
  elif [ "$iso" -eq "Zr92" ] || [ "$iso" -eq "Zr94" ]
  then
    maxrand=200
  elif [ "$iso" -eq "Zr92" ]
  then
    maxrand=250
  else
    maxrand=299
  fi
  for (( irand=0; irand<= $maxrand; irand++ ))
  do
    # Prepare the files to be modified
    cp $CurrentDir/../data/Tihange.x2m    Tihange_${iso}_${irand}.x2m
    cp $CurrentDir/../data/Tihange.access Tihange_${iso}_${irand}.access
    cp $CurrentDir/../data/Tihange.save   Tihange_${iso}_${irand}.save
    # Change the isotope from JEFF-3.1.1 to a random one from TENDL
    sed -i 's/'$iso'.*/'$iso' 1 = '$irand' ;/' Tihange_${iso}_${irand}.x2m
    # Writing Slurm instruction file
    echo "#!/bin/bash"                                                     > ${iso}_${irand}_Tihange.sh
    echo "#SBATCH -c1  ## Number of cores to be reserved"                 >> ${iso}_${irand}_Tihange.sh
    echo "#SBATCH -t 0-00:30:00"                                          >> ${iso}_${irand}_Tihange.sh
    echo "#SBATCH --mem-per-cpu=4G"                                       >> ${iso}_${irand}_Tihange.sh
    echo "srun $CurrentDir/../runVersion5.sh Tihange_${iso}_${irand}.x2m" >> ${iso}_${irand}_Tihange.sh
    # Executing the instruction file
    chmod 755 ${iso}_${irand}_Tihange.sh
    sbatch ${iso}_${irand}_Tihange.sh
  done
done

mkdir $CurrentDir/../Linux_x86_64
mv ../Linux_x86_64/* $CurrentDir/../Linux_x86_64/.

#chmod -R 777 $LaunchDir
#/bin/rm -rf $LaunchDir
