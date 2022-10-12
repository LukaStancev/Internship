#!/bin/sh
# This script has been conceived as a step-by-step guide, to help
# reproducing this work (totally or partially). Because of its high
# computation cost, it has never been executed as such, but rather
# command-by-command. So this script is perfectly untested. If executed,
# it may crash at some point (untested, again), so we...
exit
# to avoid any misuse.
# Author : V. Salino (IRSN), 10/2022

#---
#  Source environment variables
#---
. export.sh

#---
#  Download ENDF files
#---
cd ENDF
./DownloadENDF.sh
cd ..

#---
#  Produce Draglib and ACE files for Dragon and Serpent, respectively
#---
cd PyNjoy2016
cmake
make # compiles NJOY2016
cd python
./Compile_SANDY.sh
cd slurm
srun jeff3p3_shem295.sh # Best-estimate nuclear data
srun TENDL.sh # Random nuclear data, from TENDL-2019
srun SANDY.sh # Random nuclear data, sampled on-the-fly with SANDY
cd ..
./xsdir_xsdata.sh # Transform XSDIR (for MCNP, originally) in XSDATA (for Serpent)
cd ../..

#---
#  Drakkar (Dragon & Donjon) computations, with best-estimate nuclear data
#---
cd Drakkar/Version5/PyGan/src
make # compiles successively Dragon, Donjon, PyGan...
cd ../../../draglib
./convert_draglib.sh # transforms Draglib from ASCII format to binary format
cd ../data
../runVersion5.sh K-infinity.x2m
../runVersion5.sh Tihange.x2m
../runVersion5.sh Almaraz.x2m
../runVersion5.sh Fessenheim-Bugey.x2m
../runVersion5.sh Fessenheim-Bugey_DetSPH.x2m
# Exporting toward Serpent
../runVersion5.sh ExportTihToSerpent.x2m
../runVersion5.sh ExportAlmToSerpent.x2m
../runVersion5.sh ExportFshBugToSerpent.x2m

#---
#  Drakkar (Dragon & Donjon) computations, with uncertain nuclear data (randomly sampled)
#---
cd ../Slurm
srun TMC.sh
srun FSH-BUG_TMC.sh
srun FSH-BUG_TMC_with_interactions.sh
cd ../..

#---
#  D2S : feeding Serpent from Dragon datasets
#---
cd Serpent/D2S
python Assemblies_kinf.py
python FullCoreTihange.py
python FullCoreAlmaraz.py
python FullCoreFessenheimBugey.py

#---
#  Serpent computations, with best estimate nuclear data
#---
cd ../Assemblies_kinf
./SLURM.sh
cd ../FullCore # on Tihange-1 cases
./SLURM.sh
cd ../BatchHist # on Tihange-1 cases
./SLURM.sh
cd ../ParallelSpeedupBenchmark # on Tihange-1 cases
./FullCore.sh
cd ../SuperBatch
./SLURM.sh # Bugey-2, Fessenheim-1 and 2
./SLURM_Almaraz.sh

#---
#  Serpent computations, with uncertain nuclear data (randomly sampled)
#---
cd ../csw
./SLURM.sh # producing csw file to accelerate convergence for TMC uses
cd ../TMC
./GenerateInputs.sh
srun TihangeFullCore*.sh # highly intensive computations
cd ../..

#---
#  Produce plots
#---
cd Plots
./RunAll.sh
cd ..
