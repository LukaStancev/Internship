#!/bin/sh
# Usage : . export.sh

# System-specific paths (shall be modified to match with the installation of your system)
export FORTRANPATH=/soft/gcc/8.2.0/lib64/
module load anaconda3 # Load an up-to-date Python3 distribution (includes NumPy)

# Generic paths (should not require any modification)
System=`uname -s`
Sysx="`echo $System | cut -b -6`"
if [ $Sysx = "CYGWIN" ]; then
   MACH=`uname -o`
elif [ $Sysx = "AIX" ]; then
   MACH=`uname -s`
else
   MACH=`uname -sm | sed 's/[ ]/_/'`
fi
export PYTHONPATH=`pwd`/../Drakkar/Version5/PyGan/lib/$MACH/python
