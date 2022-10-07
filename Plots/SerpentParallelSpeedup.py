#
#  Plotting speedup gained with OpenMP parallel computing in Serpent
#  Usage : python3 SerpentParallelSpeedup.py
#  Author : V. Salino (IRSN), 03/2020
#

# Imports
from BasicFunctions import *
import serpentTools
from serpentTools.settings import rc
import numpy as np
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en

#---
#  Retrieve Serpent transport cycle times and compute speedup, for different
#  number of cores (from 1 to 40)
#---
times = {}
speedup = {}
respaths = glob.glob('../Serpent/ParallelSpeedupBenchmark/*.sss2_res.m')
for respath in respaths:
    res = serpentTools.read(respath)
    # Retrieve time in Serpent result file (in minutes)
    time = res.resdata['transportCycleTime'][0]
    # Retrieve case (ppm of Boron <=> control rod insertion) and number of
    # cores used directly in filename
    splitting = respath.replace('_', ' ').replace('.', ' ').split()
    ncore = int(splitting[-4])
    ppmBoron = ''.join([x for x in splitting[0] if x.isdigit()])
    if ppmBoron not in times:
        times[ppmBoron] = [[ncore, time]]
    else:
        times[ppmBoron].append([ncore, time])
for ppmBoron in times.keys():
    times[ppmBoron].sort() # by ncore (from 1 to 40)
    times[ppmBoron] = np.array(times[ppmBoron])
    speedup[ppmBoron] = np.array([times[ppmBoron][:,0], # Keep ncore
                                 times[ppmBoron][0,1]/ # time for 1 core
                                 times[ppmBoron][:,1]]) # time for ncore
#---
#  Produce the plot itself
#---
fig, ax = plt.subplots(figsize = set_size(aspect = 'speedup'))
labels = {'960': 'Groupes C et D insérés',
          '1084': 'Groupe D inséré',
          '1206': 'Toutes grappes extraites'}
for ppmBoron in labels.keys():
    ax.scatter(speedup[ppmBoron][0],
               speedup[ppmBoron][1],
               label = labels[ppmBoron], zorder = 2)
#---
#  Graph glitter
#---
ax.set(xlabel='Nombre de cœurs', ylabel="Facteur d'accélération")
ax.grid()
# Add ticks at x=1 and y=1
ax.set_xticks([1] + list(np.arange(1,10)*5))
ax.set_yticks([1] + list(np.arange(1,10)*5))
# Add an 'ideal' line and set limits accordingly to min/max
mincore = min(speedup[ppmBoron][0])
maxcore = max(speedup[ppmBoron][0])
ax.plot([mincore, maxcore], [mincore, maxcore], color = 'black',
        linestyle='--', zorder = -1)
ax.set_xlim(mincore-0.5, maxcore+0.5)
ax.set_ylim(mincore-0.5, maxcore+0.5)
# Set a non-transparent legend
plt.legend(framealpha = 1.0)
#---
#  Save plot as pdf (vectorized)
#---
fig.savefig('SerpentParallelSpeedup.pdf', bbox_inches = 'tight')
print("Plotting completed")
