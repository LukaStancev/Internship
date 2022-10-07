#!/bin/bash

for script in BatchHist.py DetectorsDiff.py Draglib.py FigurativeDistr.py Kinf.py Pearson.py PowDifferenceCP0.py PowDifferenceTih.py StandardError.py SerpentParallelSpeedup.py TMC_FSH-BUG_interactions.py TMC_FSH-BUG.py TMC_Tih_comparison.py TMC_Tih.py
do
    echo "Running "$script" ..."
    time python $script
done
