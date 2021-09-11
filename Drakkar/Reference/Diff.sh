#!/bin/bash
iexit=0

if diff -q ../Reference/JEFF33_Pdist.txt  <(grep --color=none "RELATIVE POWER" -A 16 ../Output_TIH_BestEstimate/Tihange.result | grep --color=none "  1 " -A 14) > /dev/null
then
   echo "The obtained power distributions and the reference ones are identical."
else
   echo "The obtained power distributions and the reference ones are different."
   echo "The reference power distributions are:"
   cat ../Reference/JEFF33_Pdist.txt
   echo "The obtained power distributions are:"
   grep --color=none "RELATIVE POWER" -A 16 ../Output_TIH_BestEstimate/Tihange.result | grep --color=none "  1 " -A 14
   iexit=1
fi

if diff -q ../Reference/JEFF33_keff.txt  <(grep --color=none "EFFECTIVE M" ../Output_TIH_BestEstimate/Tihange.result | cut -d"=" -f2) > /dev/null
then
   echo "The obtained k-effectives and the reference ones are identical."
else
   echo "The obtained k-effectives and the reference ones are different."
   echo "The reference k-effectives are:"
   cat ../Reference/JEFF33_keff.txt
   echo "The obtained k-effectives are:"
   grep --color=none "EFFECTIVE M" ../Output_TIH_BestEstimate/Tihange.result | cut -d"=" -f2
   iexit=1
fi

if [ "$iexit" -eq 1 ]
then
   exit 1
else
   exit 0
fi

