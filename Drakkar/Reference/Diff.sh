#!/bin/bash
iexit=0

for state in ARO D CD
do
    if diff -q ../Reference/Power${state}_TIH_JEFF33.ascii ../Output_TIH_BestEstimate/_Power${state}*.ascii
    then
        echo "The power distribution and k-effective, in state ${state}, are identical between the reference and the obtained values."
    else
        echo "The obtained and reference power distribution and/or k-effective are different. The references values are:"
        cat ../Reference/Power${state}_TIH_JEFF33.ascii
        echo "The obtained values are:"
        cat ../Output_TIH_BestEstimate/_Power${state}*.ascii
        echo "The difference(s) is/are:"
        diff ../Reference/Power${state}_TIH_JEFF33.ascii ../Output_TIH_BestEstimate/_Power${state}*.ascii
        iexit=1
    fi
done

if [ "$iexit" -eq 1 ]
then
   exit 1
else
   exit 0
fi

