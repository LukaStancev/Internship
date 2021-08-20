#!/bin/bash

for script in BatchHist.py  Draglib.py  FigurativeDistr.py  Kinf.py  PowDifference.py  StandardError.py  TMC_comparison.py  TMC.py
do
    echo "Running "$script" ..."
    time python $script
done
