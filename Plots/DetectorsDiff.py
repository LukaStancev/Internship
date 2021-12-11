#
#  Plotting detector relative difference on 2D maps
#  Usage : python3 DetectorsDiff.py
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
#  Initialize the container of all combinations of detector responses
#---
sources = ['Drakkar_DetSPH', 'Drakkar', 'Drakkar_corrected', 'Fessenheim-1',
           'Fessenheim-2', 'Bugey-2', 'ExpMean']
labels = {'Fessenheim-1'      : 'FSH1',
          'Fessenheim-2'      : 'FSH2',
          'Bugey-2'           : 'BUG2',
          'ExpMean'           : 'MoyenneExp',
          'Drakkar'           : 'Drakkar',
          'Drakkar_DetSPH'    : 'Drakkar_{DetSPH}',
          'Drakkar_corrected' : r'Drakkar_{corrig\acute{e}}'}
responses = {}
# Layout of the assemblies in the core
CoreLayout, FullCoreLayout = GetCoreLayout(157)
#---
#  Retrieve responses measured experimentally by 50 mobile fission chambers
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
    respon = pd.DataFrame(lines).replace(r'^\s*$', np.nan, regex=True)
    # Transform values from strings to floats
    respon = respon.astype('float')
    # Check that normalization is properly done, up to a criterion (0.05%)
    nelem = np.sum(respon.count())
    totsum = respon.sum().sum()
    if np.abs(totsum/nelem - 1) > 0.05/100:
        raise Exception('The normalization of the file ' + file
                        + ' is incorrect. The sum of the detector responses is'
                        + ' ' + str(totsum) + ' instead of ' + str(nelem)
                        + '. Review its detector responses.')
    # Renormalize to remove any possible small remaining residue
    respon = respon*nelem/totsum
    # Identify the source
    reactor = re.split('[/.]', file)[-2]
    responses[reactor] = np.array(respon)
#---
#  Plot relative standard deviation for each detector
#---
respon = np.array(list(responses.values()))
relstd = np.nanstd(respon, axis=0, ddof=1)/np.nanmean(respon, axis=0)*100
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
os.system('mkdir -p output_DetectorsDiff')
fig.savefig('output_DetectorsDiff/ExpRelStd.pdf', bbox_inches='tight')
# Clean-up for next plot
plt.close('all')
del fig
del ax
#---
#  Include the average of the experimental measurements
#---
totsum = 0
for reactor in responses:
    totsum = totsum + responses[reactor]
responses['ExpMean'] = totsum/len(responses)
#---
#  Drakkar detector responses
#---
# Retrieve and deploy in 2D the Drakkar responses
link = '../Drakkar/Output_FSH-BUG_BestEstimate/_PowerARO_*.ascii'
respon = Deploy2D(loadlcm(link)['RESPON'], FullCoreLayout)
# Filter out unmeasured responses with NaNs, in order to have a
# normalization identical to the measurements, i.e.
# sum(50 measured positions)=50.
respon_filtered = np.where(np.isnan(responses['Fessenheim-2']), np.nan, respon)
normfactor = np.nanmean(respon_filtered)
# Filter out the reflectors (=0 in FullCoreLayout) with NaNs and apply that
# normalization
responses['Drakkar'] = np.where(FullCoreLayout == 0, np.nan, respon)/normfactor
#---
#  Drakkar detector responses, computed from SPH-corrected detector cross
#  sections (not recommended ; retrieved here for demonstration purposes)
#---
# Retrieve and deploy in 2D the Drakkar responses
link = '../Drakkar/Output_FSH-BUG_BestEstimate_DetSPH/_PowerARO_*.ascii'
respon = Deploy2D(loadlcm(link)['RESPON'], FullCoreLayout)
# Filter out unmeasured responses with NaNs, in order to have a
# normalization identical to the measurements, i.e.
# sum(50 measured positions)=50.
respon_filtered = np.where(np.isnan(responses['Fessenheim-2']), np.nan, respon)
normfactor = np.nanmean(respon_filtered)
# Filter out the reflectors (=0 in FullCoreLayout) with NaNs and apply that
# normalization
responses['Drakkar_DetSPH'] = (np.where(FullCoreLayout == 0, np.nan, respon)/
                               normfactor)
#---
#  Correct the deterministic bias of Drakkar's detector responses, based on the
#  differences in power distribution between Drakkar and Serpent.
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
responses['Drakkar_corrected'] = (responses['Drakkar']*powerSerpent
                                                      /powerDrakkar)
#---
#  Plot all the responses
#---
os.system('mkdir -p output_DetectorsDiff')
amin = math.inf
amax = -math.inf
for source in sources:
    amin = min(amin, np.nanmin(responses[source]))
    amax = max(amax, np.nanmax(responses[source]))
for source in sources:
    respon = responses[source]
    # Produce the plot
    fig, ax = plt.subplots(1, figsize = set_size('square'))
    axins = inset_axes(ax, width = "100%", height = "3%",
                       loc = 'lower center',
                       bbox_to_anchor = (0.05, -0.17, 0.95, 1.0),
                       bbox_transform = ax.transAxes, borderpad = 5)
    im = ax.imshow(respon, interpolation = 'nearest', cmap = 'coolwarm',
                   vmin = amin, vmax = amax)
    cbar = plt.colorbar(im, cax = axins, orientation = 'horizontal')
    current_cmap = im.get_cmap()
    current_cmap.set_bad(color = 'white')
    cbar.ax.tick_params(labelsize = 12)
    VisualAids(ax)
    # In each assembly, add the axially-integrated reponses
    for i in range(len(respon)):
        for j in range(len(respon)):
            if not np.isnan(respon[i, j]):
                text = str('{:.3f}'.format(round(respon[i, j], 3)))
                ax.text(j, i, text, color = 'black',
                        ha = 'center', va = 'center')
    # Save plot as pdf (vectorized)
    FullCoreGlitter(ax, FullCoreLayout)
    fig.savefig('output_DetectorsDiff/' + source + '.pdf',
                bbox_inches='tight')
    # Clean-up for next plot
    plt.close('all')
    del fig
    del ax
#---
#  Compute relative differences between all sources and plot them
#---
for source in sources:
    responses[source] = pd.DataFrame(responses[source])
# Get every two-by-two combinations between every sources
for (a, b) in list(itertools.combinations(sources, 2)):
    discrp = (responses[a]/responses[b] - 1)*100
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
    os.system('mkdir -p output_DetectorsDiff')
    fig.savefig('output_DetectorsDiff/' + a + '_vs_' + b + '.pdf',
                bbox_inches='tight')
    # Clean-up for next plot
    plt.close('all')
    del fig
    del ax
