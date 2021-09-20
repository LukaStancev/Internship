#
#  Plotting detector relative difference on 2D maps
#  Usage : python3 ActivitiesDiff.py
#  Author ; V. Salino (IRSN), 09/2021
#

# Imports
import lcm
import numpy as np
import pandas as pd
from scipy import stats
import os
import math
import subprocess
import glob
import re
import csv
import string
import itertools
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from BasicFunctions import *
import serpentTools
from serpentTools.settings import rc
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en

#---
#  Initialize the container of all combinations of detector activities
#---
activities = {}
# Layout of the assemblies in the core
CoreLayout, FullCoreLayout = GetCoreLayout(157)
FullCoreLayout4X = np.repeat(np.repeat(FullCoreLayout,
                                       2, axis = 0), 2, axis = 1)
#---
#  Retrieve activities measured experimentally by 50 mobile fission chambers
#---
for file in (set(glob.glob('../Measurements/CP0/*'))
# Minus Almaraz, which contains only measurements-derived power distribution
             - set(glob.glob('../Measurements/CP0/Almaraz*'))):
    print('Reading ' + file)
    lines = []
    with open(file) as csvfile:
        reader = csv.reader(decomment(csvfile))
        for line in reader:
            lines.append(line)
    # Replace spaces with NaNs
    activ = pd.DataFrame(lines).replace(r'^\s*$', np.nan, regex=True)
    # Transform values from strings to floats
    activ = activ.astype('float')
    # Check that normalization is properly done, up to a criterion (0.05%)
    nelem = np.sum(activ.count())
    totsum = activ.sum().sum()
    if np.abs(totsum/nelem - 1) > 0.05/100:
        raise Exception('The normalization of the file ' + file
                        + ' is incorrect. The sum of the activities is '
                        + str(totsum) + ' instead of ' + str(nelem)
                        + '. Review its activities.')
    # Renormalize to remove any possible small remaining residue
    activ = activ*nelem/totsum
    # Identify the source
    reactor = re.split('[/.]', file)[-2]
    activities[reactor] = np.array(activ)
#---
#  Include the average of the experimental measurements
#---
totsum = 0
for reactor in activities:
    totsum = totsum + activities[reactor]
activities['ExpMean'] = totsum/len(activities)
#---
#  Retrieve best-estimate Drakkar and Serpent power distributions and also
#  CASMO-5 detector-to-power ratio
#---
path = '../Drakkar/Output_FSH-BUG_BestEstimate/'
link = path + '_PowerARO_4608cfa17412e538f9d0c24af965a2f4.ascii'
power = loadlcm(link)['POWER-CHAN']
# Each FullCoreLayout value is duplicated in 4 values
power = matrixaverage(Deploy2D(power, FullCoreLayout4X), 2)
# Filter out reflectors (=0 in FullCoreLayout) with NaNs
power = np.where(FullCoreLayout == 0, np.nan, power)
# Normalize
powerDrakkar = power/np.nanmean(power)
# Retrieve Serpent power distribution
links = glob.glob('../Serpent/SuperBatch/Fessenheim-Bugey*.sss2_det0.m')
# Storage for independant runs' results
indruns = []
for link in links:
    indruns.append(ReadPdistrSerpent(link, FullCoreLayout))
powerSerpent = Deploy2D(np.mean(indruns, axis = 0), FullCoreLayout)
DetToPow = np.array(pd.read_csv('CASMO5_DetectorToPowerRatio.csv', sep = ";",
                    header=None))
DetToPow = UnfoldQuarter(np.rot90(DetToPow))
#---
#  Retrieve randomly sampled Drakkar powers and activities
#---
activitiesTMC = []
# Short execution:
#for link in glob.glob("../Drakkar/Output_FSH-BUG_TMC/_PowerARO_1*"):
for link in glob.glob("../Drakkar/Output_FSH-BUG_TMC/_PowerARO_*"):
    print(link)
    # Retrieve power
    power = matrixaverage(Deploy2D(loadlcm(link)['POWER-CHAN'],
                                   FullCoreLayout4X), 2)
    power = np.where(FullCoreLayout == 0, np.nan, power)
#    # Retrieve and deploy in 2D the Drakkar activities
#    activ = Deploy2D(loadlcm(link)['RESPON'], FullCoreLayout)
    # Until the Drakkar detector problems are resolved, combine Drakkar power
    # distribution and CASMO-5 detector results. In addition, correct for
    # deterministic bias using Serpent.
    activ = power*DetToPow*powerSerpent/powerDrakkar
    # Filter out unmeasured activities with NaNs, in order to have a
    # normalization identical to the measurements, i.e.
    # sum(50 measured positions)=50.
    activ_filtered = np.where(np.isnan(activities['Fessenheim-2']), np.nan,
                              activ)
    normfactor = np.nanmean(activ_filtered)
    # Filter out the reflectors (=0 in FullCoreLayout) with NaNs and apply that
    # normalization
    activ = np.where(FullCoreLayout == 0, np.nan, activ)/normfactor
    activitiesTMC.append(activ)
# Convert this list of 1D numpy array into a 2D numpy array
activitiesTMC = np.array(activitiesTMC)
#---
#  Statistics
#---
ndet = np.count_nonzero(~np.isnan(activities['ExpMean']))
uncertainty = 0.02 # 2%
ncalc = len(activitiesTMC)
khisqrs = np.empty(ncalc)
i = 0
for i in range(0, ncalc):
    khisqrs[i] = np.nansum(((activitiesTMC[i] - activities['ExpMean'])
                             /uncertainty)**2)/ndet
print(stats.describe(khisqrs))
# BMC and BFMC weights :
# * Rochman, D., Bauge, E., Vasiliev, A. et al. "Monte Carlo nuclear data
#   adjustment via integral information." Eur. Phys. J. Plus 133, 537 (2018).
#   https://tendl.web.psi.ch/bib_rochman/correlation.pdf
# * Rochman, D., Vasiliev, A., Ferroukhi, H. et al. "Correlation nu_p-sigma for
#   U-Pu in the thermal and resonance neutron range via integral information."
#   Eur. Phys. J. Plus 134, 453 (2019).
#   https://tendl.web.psi.ch/bib_rochman/epjp1900495-offprints.pdf
weights_BMC = np.empty_like(khisqrs)
weights_BFMC = np.empty_like(khisqrs)
weights_BMC = np.exp(-khisqrs/2)
weights_BFMC = np.exp(-(khisqrs/np.min(khisqrs))**2)
print(stats.describe(weights_BMC))
print(stats.describe(weights_BFMC))
#---
#  Plotting
#---
# Initialize figure
fig, axs = plt.subplots(len(FullCoreLayout), len(FullCoreLayout),
                        sharex = 'all',
                        figsize = set_size('square'),
                        gridspec_kw = {'hspace': 0, 'wspace': 0})
for x in range(0, len(FullCoreLayout[0, :])):
    for y in range(0, len(FullCoreLayout[:, 0])):
        if FullCoreLayout[x, y] == 0:
            axs[x, y].set_axis_off()
        else:
            # Plot an histogram with the assembly power
            axs[x, y].hist(activitiesTMC[:, x, y], bins = 20,
                           density = True)
            # Compute relative standard deviation for each assembly
            mean = np.mean(activitiesTMC[:, x, y])
            relstd = np.std(activitiesTMC[:, x, y], ddof = 1)/mean*100
            textstr = (r'$\mu$=$%.2f$' % mean + '\n'
                       + r'$\sigma$=$%.1f$' % (relstd, ) + '\%')
            # Print statistics for each assembly
            axs[x, y].text(0.05, 0.95, textstr, fontsize = 6,
                           transform=axs[x, y].transAxes,
                           verticalalignment = 'top')
            # Remove Y-ticks (unnecessary for a probability density)
            axs[x, y].set_yticks([])
for x in range(0, len(CoreLayout[0, :])):
    for y in range(0, len(CoreLayout[:, 0])):
        # Assemblies are squares so we'd like square subplots
        axs[x, y].set_aspect(1./axs[x, y].get_data_ratio())
        # Remove some of the axes, in order to avoid doubled axes
        if y != len(CoreLayout[:, 0]) - 1:
            # If there is a graph to its right, do not draw the right
            # axis
            if CoreLayout[x, y + 1] == 1:
                axs[x, y].spines['right'].set_visible(False)
        if x != len(CoreLayout[:, 0]) - 1:
            # If there is a graph at his bottom, do not draw the bottom
            # axis
            if CoreLayout[x + 1, y] == 1:
                axs[x, y].spines['bottom'].set_visible(False)
# Add a title and save plot as pdf (vectorized)
os.system('mkdir -p output_TMC_FSH-BUG')
fig.savefig('output_TMC_FSH-BUG/TMC.pdf', bbox_inches='tight')
# Clean-up for next plot
plt.close('all')
