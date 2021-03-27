#!/bin/bash
# Download in the execution directory all the ENDF files used within this project
# Usage : ./DownloadENDF.sh
#---
#  JEFF-3.1.1
#---
mkdir -p JEFF-3.1.1
cd JEFF-3.1.1
# Neutron data
wget https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_31/JEFF311/FINAL/JEFF311N_0_IND.zip
# Thermal scattering law (TSL) data
wget https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_31/JEFF31/JEFF31TS_INDIV.tar.gz
# Radioactive decay data (RDD) and neutron-induced fission yields (NFY)
wget https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_31/JEFF-311_RDD/JEFF311RDD_ALL.OUT
wget https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_31/JEFF311/JEFF-311FY/JEFF311NFY.ASC
# Extract data
unzip JEFF311N_0_IND.zip -d JEFF311N_0_IND
mkdir -p JEFF31TS_INDIV
tar -xvzf JEFF31TS_INDIV.tar.gz -C JEFF31TS_INDIV
rm -f JEFF311N_0_IND.zip JEFF31TS_INDIV.tar.gz
cd ..
#---
#  JEFF-3.1.2
#---
mkdir -p JEFF-3.1.2
cd JEFF-3.1.2
# Neutron data
wget https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_31/JEFF312/final/JEFF312N.zip
# Extract data
unzip JEFF312N.zip -d JEFF312N
rm -f JEFF312N.zip
cd ..
#---
#  JEFF-3.2
#---
mkdir -p JEFF-3.2
cd JEFF-3.2
# Neutron data
wget https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_32/JEFF32N.tar.gz
# Extract data
mkdir -p JEFF32N
tar -xvzf JEFF32N.tar.gz -C JEFF32N
rm -f JEFF32N.tar.gz
cd ..
#---
#  JEFF-3.3
#---
mkdir -p JEFF-3.3
cd JEFF-3.3
# Neutron data
wget https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/JEFF33-n.tgz
# Thermal scattering law (TSL) data
wget https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/JEFF33-tsl.tgz
# Radioactive decay data (RDD) and neutron-induced fission yields (NFY)
wget https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/JEFF33-rdd_all.asc
wget https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/JEFF33-nfy.asc
# Extract data
mkdir -p JEFF33-n JEFF33-tsl
tar -xvzf JEFF33-n.tgz -C JEFF33-n
tar -xvzf JEFF33-tsl.tgz
rm -f JEFF33-n.tgz JEFF33-tsl.tgz
mv JEFF33-n/endf6/* JEFF33-n/.
rm -rf JEFF33-n/endf6
cd ..
#---
#  TENDL-2019
#---
mkdir -p TENDL-2019
cd TENDL-2019
wget https://tendl.web.psi.ch/tendl_2019/tar_files/Fe056.random.tgz
mkdir -p Fe56
tar -xvzf Fe056.random.tgz -C Fe56
rm -f Fe056.random.tgz
gzip -d Fe56/*.gz
for iso in Zr090 Zr091 Zr092 Zr094 Zr096 Ag107 Ag109 In115 Cd106 Cd108 Cd110 Cd111 Cd112 Cd113 Cd114 Cd116 Fe054 Cr052 Ni058
do
    # Strip the numbers to keep only the name of the element (Zr, Ag, ...)
    element=$(echo $iso | sed 's/[0-9]//g')
    Z=$(echo $iso | sed 's/[A-Za-z]//g')
    # Strip leading zero (Zr090 -> Zr90)
    Z_no0=$(echo $Z | sed 's/^0//')
    iso_no0="$element$Z_no0"
    mkdir -p $iso_no0
    # Download
    wget https://tendl.web.psi.ch/tendl_2019/neutron_file/$element/$iso/lib/endf/random/random.$iso.tgz --directory-prefix=./$iso_no0
    # Extract
    tar -xvzf $iso_no0/random.$iso.tgz -C $iso_no0
    rm -f $iso_no0/random.$iso.tgz
    # Repatriate files that are scattered in a swarm of directories
    find $iso_no0/ -type f | xargs -I '{}' mv '{}' $iso_no0/.
    rm -rf $iso_no0/$iso
    # Uncompress
    gzip -d $iso_no0/*.gz
done
cd ..
