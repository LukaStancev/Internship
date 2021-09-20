#
#  Plotting detector relative difference on 2D maps
#  Usage : python3 ActivitiesDiff.py
#  Author ; V. Salino (IRSN), 09/2021
#

# Imports
import lcm
import numpy as np
import pandas as pd
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
sources = ['DrakC5', 'SerpC5', 'Drakkar', 'Drakkar_corrected', 'Fessenheim-1', 'Fessenheim-2',
           'Bugey-2', 'ExpMean']
labels = {'Fessenheim-1'      : 'FSH1',
          'Fessenheim-2'      : 'FSH2',
          'Bugey-2'           : 'BUG2',
          'ExpMean'           : 'MoyenneExp',
          'Drakkar'           : 'Drakkar',
          'DrakC5'            : 'DrakC5',
          'SerpC5'            : 'SerpC5',
          'Drakkar_corrected' : r'Drakkar_{corrig\acute{e}}'}
activities = {}
# Layout of the assemblies in the core
CoreLayout, FullCoreLayout = GetCoreLayout(157)
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
#  Plot relative standard deviation for each detector
#---
activ = np.array(list(activities.values()))
relstd = np.nanstd(activ, axis=0, ddof=1)/np.nanmean(activ, axis=0)*100
# Produce the plot
fig, ax = plt.subplots(1, figsize = set_size('square'))
im = ax.imshow(relstd, interpolation = 'nearest', cmap = 'coolwarm',
               vmin = 0.3, vmax = 2.8)
# Add a colorbar below
axins = inset_axes(ax, width = "100%", height = "3%",
                   loc = 'lower center',
                   bbox_to_anchor = (0.05, -0.17, 0.95, 1.0),
                   bbox_transform = ax.transAxes, borderpad = 5)
cbar = plt.colorbar(im, cax = axins, orientation = 'horizontal',
                    ticks = list(0.5 + np.arange(0, 5)/2))
current_cmap = im.get_cmap()
current_cmap.set_bad(color = 'white')
# Add labels to the colorbar
cbar.ax.tick_params(labelsize = 12)
cbarlabels = [str(i) + '\%' for i in list(0.5 + np.arange(0, 5)/2)]
cbar.ax.set_xticklabels(cbarlabels)
VisualAids(ax)
# Add relstd inside each box
nprelstd = np.array(relstd)
for i in range(len(nprelstd)):
    for j in range(len(nprelstd)):
        if not np.isnan(nprelstd[i, j]):
            text = str('{:.1f}'.format(round(relstd[i, j], 2)))
            ax.text(j, i, text, color = 'black',
                    ha = 'center', va = 'center')
# Save plot as pdf (vectorized)
FullCoreGlitter(ax, FullCoreLayout)
os.system('mkdir -p output_ActivitiesDiff')
fig.savefig('output_ActivitiesDiff/ExpRelStd.pdf', bbox_inches='tight')
# Clean-up for next plot
plt.close('all')
del fig
del ax
#---
#  Include the average of the experimental measurements
#---
totsum = 0
for reactor in activities:
    totsum = totsum + activities[reactor]
activities['ExpMean'] = totsum/len(activities)
#---
#  Drakkar detector activities
#---
# Retrieve and deploy in 2D the Drakkar activities
link = glob.glob('../Drakkar/Output_FSH-BUG_BestEstimate/_PowerARO*')[0]
activ = Deploy2D(loadlcm(link)['RESPON'], FullCoreLayout)
# Filter out unmeasured activities with NaNs, in order to have a
# normalization identical to the measurements, i.e.
# sum(50 measured positions)=50.
activ_filtered = np.where(np.isnan(activities['Fessenheim-2']), np.nan, activ)
normfactor = np.nanmean(activ_filtered)
# Filter out the reflectors (=0 in FullCoreLayout) with NaNs and apply that
# normalization
activities['Drakkar'] = np.where(FullCoreLayout == 0, np.nan, activ)/normfactor
#---
#  Correct the deterministic bias of Drakkar's detector activities, based on
#  the differences in power distribution between Drakkar and Serpent.
#---
# Retrive Drakkar power distribution
power = loadlcm(link)['POWER-CHAN']
# Each FullCoreLayout value is duplicated in 4 values
FullCoreLayoutX4 = np.repeat(np.repeat(FullCoreLayout, 2, axis=0), 2, axis=1)
power = matrixaverage(Deploy2D(power, FullCoreLayoutX4), 2)
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
activities['Drakkar_corrected'] = (activities['Drakkar']*powerSerpent
                                                        /powerDrakkar)
#---
#  CASMO-5
#---
DetToPow = np.array(pd.read_csv('CASMO5_DetectorToPowerRatio.csv', sep = ";",
                    header=None))
DetToPow = UnfoldQuarter(np.rot90(DetToPow))
# Apply these factors to Drakkar and Serpent
for [source, power] in [['DrakC5', powerDrakkar], ['SerpC5', powerSerpent]]:
    activ = power*DetToPow
    # Filter out unmeasured activities with NaNs, in order to have a
    # normalization identical to the measurements, i.e.
    # sum(50 measured positions)=50.
    activ_filtered = np.where(np.isnan(activities['Fessenheim-2']), np.nan, activ)
    normfactor = np.nanmean(activ_filtered)
    # Filter out the reflectors (=0 in FullCoreLayout) with NaNs and apply that
    # normalization
    activities[source] = activ/normfactor
#---
#  Plot all the activities
#---
os.system('mkdir -p output_ActivitiesDiff')
amin = math.inf
amax = -math.inf
for source in sources:
    amin = min(amin, np.nanmin(activities[source]))
    amax = max(amax, np.nanmax(activities[source]))
for source in sources:
    activ = activities[source]
    # Produce the plot
    fig, ax = plt.subplots(1, figsize = set_size('square'))
    axins = inset_axes(ax, width = "100%", height = "3%",
                       loc = 'lower center',
                       bbox_to_anchor = (0.05, -0.17, 0.95, 1.0),
                       bbox_transform = ax.transAxes, borderpad = 5)
    im = ax.imshow(activ, interpolation = 'nearest', cmap = 'coolwarm',
                   vmin = amin, vmax = amax)
    cbar = plt.colorbar(im, cax = axins, orientation = 'horizontal')
    current_cmap = im.get_cmap()
    current_cmap.set_bad(color = 'white')
    cbar.ax.tick_params(labelsize = 12)
    VisualAids(ax)
    # In each assembly, add the axially-integrated activity
    for i in range(len(activ)):
        for j in range(len(activ)):
            if not np.isnan(activ[i, j]):
                text = str('{:.3f}'.format(round(activ[i, j], 3)))
                ax.text(j, i, text, color = 'black',
                        ha = 'center', va = 'center')
    # Save plot as pdf (vectorized)
    FullCoreGlitter(ax, FullCoreLayout)
    fig.savefig('output_ActivitiesDiff/' + source + '.pdf',
                bbox_inches='tight')
    # Clean-up for next plot
    plt.close('all')
    del fig
    del ax
#---
#  Compute relative differences between all sources and plot them
#---
for source in sources:
    activities[source] = pd.DataFrame(activities[source])
# Get every two-by-two combinations between every sources
for (a, b) in list(itertools.combinations(sources, 2)):
    discrp = (activities[a]/activities[b] - 1)*100
    print(discrp.min().min(), discrp.max().max())
    # Produce the plot
    fig, ax = plt.subplots(1, figsize = set_size('square'))
    im = ax.imshow(discrp, interpolation = 'nearest', cmap = 'coolwarm',
                   vmin = -10, vmax = 10)
    # Add a colorbar below
    axins = inset_axes(ax,
                       width = "100%",
                       height = "3%",
                       loc = 'lower center',
                       bbox_to_anchor = (0.05, -0.17, 0.95, 1.0),
                       bbox_transform = ax.transAxes,
                       borderpad = 5
                       )
    cbar = plt.colorbar(im, cax = axins, orientation = 'horizontal',
                        ticks = list(np.arange(-2, 3)*5))
    current_cmap = im.get_cmap()
    current_cmap.set_bad(color = 'white')
    # Add labels to the colorbar
    cbar.ax.tick_params(labelsize = 12)
    cbarlabels = [str(i) + '\%' for i in list(np.arange(-2, 3)*5)]
    cbar.ax.set_xticklabels(cbarlabels)
    VisualAids(ax)
    # Add percent deviation inside each box
    npdiscrp = np.array(discrp)
    for i in range(len(npdiscrp)):
        for j in range(len(npdiscrp)):
            if not np.isnan(npdiscrp[i, j]):
                text = str('{:.1f}'.format(round(npdiscrp[i, j], 2)))
                ax.text(j, i, text, color = 'black',
                        ha = 'center', va = 'center')
    # Add formula
    text = (r'$\frac{\mathrm{' + labels[a] + '}-\mathrm{' + labels[b] + '}}'
            + '{\mathrm{' + labels[b] + '}} (\%)$')
    ax.text(2.5, 14.2, text, fontsize = 12, color = 'black',
            ha = 'center', va = 'center')
    # Add a title
    if lang == 'en':
        fig.suptitle('Relative difference in detector between\n'
                     + a + ' and ' + b + ' on \nfirst zero-power '
                     'start-up')
    # Save plot as pdf (vectorized)
    FullCoreGlitter(ax, FullCoreLayout)
    os.system('mkdir -p output_ActivitiesDiff')
    fig.savefig('output_ActivitiesDiff/' + a + '_vs_' + b + '.pdf',
                bbox_inches='tight')
    # Clean-up for next plot
    plt.close('all')
    del fig
    del ax
