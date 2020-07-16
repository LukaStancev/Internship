#!/bin/bash
System=`uname -s`
Sysx="`echo $System | cut -b -6`"
if [ $Sysx = "CYGWIN" ]; then
   MACH=`uname -o`
elif [ $Sysx = "AIX" ]; then
   MACH=`uname -s`
else
   MACH=`uname -sm | sed 's/[ ]/_/'`
fi

iexit=0

if diff -q ../Reference/JEFF311_Pdist.txt  <(grep --color=none "RELATIVE POWER" -A 16 ../"$MACH"/JEFF311.result | grep --color=none "  1 " -A 14) > /dev/null
then
   echo "The obtained power distributions and the reference ones are identical."
else
   echo "The obtained power distributions and the reference ones are different."
   echo "The reference power distributions are:"
   cat ../Reference/JEFF311_Pdist.txt
   echo "The obtained power distributions are:"
   grep --color=none "RELATIVE POWER" -A 16 ../"$MACH"/JEFF311.result | grep --color=none "  1 " -A 14
   iexit=1
fi

if diff -q ../Reference/JEFF311_keff.txt  <(grep --color=none "EFFECTIVE M" ../"$MACH"/JEFF311.result | cut -d"=" -f2) > /dev/null
then
   echo "The obtained k-effectives and the reference ones are identical."
else
   echo "The obtained k-effectives and the reference ones are different."
   echo "The reference k-effectives are:"
   cat ../Reference/JEFF311_keff.txt
   echo "The obtained k-effectives are:"
   grep --color=none "EFFECTIVE M" ../"$MACH"/JEFF311.result | cut -d"=" -f2
   iexit=1
fi

if [ "$iexit" -eq 1 ]
then
   exit 1
else
   exit 0
fi

