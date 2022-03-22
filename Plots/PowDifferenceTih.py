#
#  Plotting power relative difference on 2D maps, for Tihange
#  Usage : python3 PowDifferenceTih.py
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import lcm
import numpy as np
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
#  This function is intended to operate with FuncFormatter, to add a permil
#  character at the end of each tick label
#---
def permil(x, pos):
    return str('{:.1f}'.format(x)) + r'\textperthousand'

#---
#  Initialize the container of all combinations of power distributions
#---
controlrods = ['ARO', 'D', 'CD']
CRtext = {}
CRtext['ARO'] = 'all rods out'
CRtext['D'] = 'D rod bank inserted'
CRtext['CD'] = 'C and D rod banks inserted'
sources = ['Drakkar', 'Serpent (fission)', 'Serpent', 'EDF', 'Framatome']
powers = {}
for controlrod in controlrods:
    powers[controlrod] = {}
#---
#  Measurement-derived power distributions
#---
for file in list(glob.glob('../Measurements/Tihange-1/*.csv')):
    meas = []
    with open(file) as csvfile:
        reader = csv.reader(decomment(csvfile))
        for row in reader:
            meas.append(np.array(row).astype(np.float))
    # Transform a list of 1D NumPy arrays of irregular lengths into a 2D NumPy
    # array
    MeasEighth = np.zeros([len(meas), len(max(meas, key = lambda x: len(x)))])
    for i,j in enumerate(meas):
        MeasEighth[i][0:len(j)] = j
    # Unfold the eighth map into a full map
    meas = UnfoldEighth(MeasEighth)
    # Check that normalization is properly done, up to a criterion (0.1%)
    if np.abs(np.sum(meas)/np.count_nonzero(meas) - 1) > 0.001:
        raise Exception('The normalization of the file ' + file
                        + ' is incorrect. The sum of the powers is '
                        + str(np.sum(meas)) + ' instead of '
                        + str(np.count_nonzero(meas)) + '. Review its powers.')
    # Remove all powers equal to zero
    meas = meas[np.nonzero(meas)]
    # Identify the source and the control rod insertion
    basename = file.split('/')[-1].split('.')[0]
    if len(basename.split('_')) != 2:
        raise Exception(file + ' unexpected.')
    basename = basename.split('_')
    controlrod = basename[0]
    if basename[1] == 'EDF':
        powers[controlrod]['EDF'] = meas
    elif basename[1] == 'Framatome':
        powers[controlrod]['Framatome'] = meas
    else:
        raise Exception(basename[1] + ' unknown.')
# Layout of the assemblies in the core
CoreLayout, FullCoreLayout = GetCoreLayout(len(powers['ARO']['EDF']))
#---
#  Serpent power distributions
#---
serpentfiles = glob.glob('../Serpent/BatchHist/Tihange*ppm_*.sss2_his0.m')
relstd = {}
for reaction in ['-8', '-80']:
    # Storage for independant runs' results
    indruns = {}
    for controlrod in controlrods:
        indruns[controlrod] = []
    for file in list(serpentfiles):
        cbor = file.split('/')[-1].split('_')[0]
        cbor = int(re.findall(r'\d+', cbor)[0])
        if cbor == 1206:
            controlrod = 'ARO'
        elif cbor == 1084:
            controlrod = 'D'
        elif cbor == 960:
            controlrod = 'CD'
        else:
            raise Exception(cbor + ' unknown.')
        # Read that Serpent history file
        indruns[controlrod].append(ReadPdistrHistSerpent(file, FullCoreLayout,
                                                         reaction = reaction))
    for controlrod in controlrods:
        # Convert this list of n-D NumPy array into a (n+1)-D NumPy array
        indruns[controlrod] = np.array(indruns[controlrod])
        # Mean indruns[runs, generations, iassembly] over active generations
        indact = np.mean(indruns[controlrod][:, 198: ,:], axis = 1)
        # Mean over independant runs
        mean = np.mean(indact, axis = 0)
        if reaction == '-80':
            powers[controlrod]['Serpent'] = mean
            # Relative standard error of the mean (‰)
            std = np.std(indact, axis = 0, ddof = 1)
            relstd[controlrod] = std/np.sqrt(np.shape(indact)[0])/mean*1000
        elif reaction == '-8':
            powers[controlrod]['Serpent (fission)'] = mean
#---
#  Drakkar power distributions
#---
for controlrod in controlrods:
    os.system('ln -s ../Drakkar/Output_TIH_BestEstimate/_Power' + controlrod
              + '_*.ascii .')
    for file in list(glob.glob('_Power' + controlrod + '_*.ascii')):
        # Loading data while removing beginning-of-filename underscore
        ResultFile = lcm.new('LCM_INP', file[1:])
        os.system("rm " + file)
        power = ResultFile['POWER-CHAN']
        # Compute mean in order to normalize power to 1.0 per average assembly
        mean = np.sum(power)/len(power)
        power = power/mean
        powers[controlrod]['Drakkar'] = power
#---
#  Plot all the power distributions and also Serpent standard errors
#---
os.system('mkdir -p output_PowDifferenceTih')
pmin = math.inf
pmax = -math.inf
relstdmin = math.inf
relstdmax = -math.inf
for controlrod in controlrods:
    for source in sources:
        pmin = min(pmin, np.min(powers[controlrod][source]))
        pmax = max(pmax, np.max(powers[controlrod][source]))
    if 'Serpent' in sources:
        relstdmin = min(relstdmin, np.min(relstd[controlrod]))
        relstdmax = max(relstdmax, np.max(relstd[controlrod]))
        sourcesinnextloop = sources + ['Serpent_SE']
    else:
        sourcesinnextloop = sources
for controlrod in controlrods:
    for source in sourcesinnextloop:
        # Show only a quarter of the core
        if source == 'Serpent_SE':
            plot = fold(Deploy2D(relstd[controlrod], FullCoreLayout))
            minplot = relstdmin
            maxplot = relstdmax
        else:
            plot = fold(Deploy2D(powers[controlrod][source], FullCoreLayout))
            minplot = pmin
            maxplot = pmax
        # Replace every zeroes with NaNs, in order to easily get them colored
        # in white (with set_bad)
        plot = np.where(np.isclose(plot, 0), float('nan'), plot)
        # Produce the plot
        fig, ax = plt.subplots(1, figsize = set_size('halfsquare'))
        axins = inset_axes(ax,
                           width = "3%",
                           height = "100%",
                           loc = 3,
                           bbox_to_anchor = (1.02, 0., 1, 1),
                           bbox_transform = ax.transAxes,
                           borderpad = 0)
        im = ax.imshow(plot, interpolation = 'nearest', cmap = 'coolwarm',
                       vmin = minplot, vmax = maxplot)
        cbar = plt.colorbar(im, cax = axins, orientation = 'vertical')
        if controlrod == 'CD':
            current_cmap = im.get_cmap()
            current_cmap.set_bad(color = 'white')
            cbar.ax.tick_params(labelsize = 12)
            if source == 'Serpent_SE':
                cbar.ax.yaxis.set_major_formatter(FuncFormatter(permil))
        else:
            cbar.remove()
        # In each assembly, add the power or its Serpent standard error
        for i in range(len(plot)):
            for j in range(len(plot)):
                if not np.isnan(plot[i, j]):
                    text = str('{:.2f}'.format(round(plot[i, j], 2)))
                    ax.text(j, i, text, color = 'black',
                            ha = 'center', va = 'center')
        # Add Battleship-style coordinates
        ax.xaxis.tick_bottom()
        xlabels = list(string.ascii_uppercase)[:8]
        xlabels.reverse()
        ylabels = list(np.arange(1, 9))
        plt.setp(ax, xticks = np.arange(8), xticklabels = xlabels,
                 yticks = np.arange(8), yticklabels = ylabels)
        ax.tick_params(axis = 'both', which = 'both', length = 0)
        #---
        #  Add a title and save plot as pdf (vectorized)
        #---
        if lang == 'en':
            if source == 'EDF' or source == 'Framatome':
                txt = source + ' pseudo-measured'
            else:
                txt = source + ' JEFF-3.3'
            txt = (txt + ' power distribution on Tihange\nfirst zero-power '
                         + 'start-up, ' + CRtext[controlrod])
            if source == 'Serpent':
                nbruns = np.shape(indruns[controlrod])[0]
                txt += ('\n' + r'Relative standard error of the mean '
                        + '($1\sigma$) from ' + str(nbruns) + '\n'
                        + 'independant simulations with different random seeds')
                plt.subplots_adjust(top = 0.85)
            st = fig.suptitle(txt)
            fig.savefig('output_PowDifferenceTih/' + source + '_' + controlrod
                        + '.pdf', bbox_inches='tight', extra_artists=[st])
        elif lang == 'fr':
            fig.savefig('output_PowDifferenceTih/' + source + '_' + controlrod
                        + '.pdf', bbox_inches='tight')
        #---
        #  Clean-up for next plot
        #---
        plt.close('all')
        del fig
        del ax
#---
#  Compute relative differences between all sources and plot them
#---
for controlrod in controlrods:
    # Get every two-by-two combinations between every sources
    for (a, b) in list(itertools.combinations(sources, 2)):
        # Compute relative deviation (discrepancy) between the two selected
        # sources ('a' and 'b')
        discrp1D = (powers[controlrod][a]/powers[controlrod][b] - 1)*100
        discrp = Deploy2D(discrp1D, FullCoreLayout)
        # Show only a quarter of the core
        discrp = fold(discrp)
        # Replace every zeroes with NaNs, in order to easily get them colored
        # in white (with set_bad)
        discrp = np.where(np.isclose(discrp, 0), float('nan'), discrp)
        # Produce the plot
        fig, ax = plt.subplots(1, figsize = set_size('halfsquare'))
        axins = inset_axes(ax,
                           width = "3%",
                           height = "100%",
                           loc = 3,
                           bbox_to_anchor = (1.02, 0., 1, 1),
                           bbox_transform = ax.transAxes,
                           borderpad = 0)
        im = ax.imshow(discrp, interpolation = 'nearest', cmap = 'coolwarm',
                       vmin = -10, vmax = 10)
        cbar = plt.colorbar(im, cax = axins, orientation = 'vertical')
        if controlrod == 'CD':
            current_cmap = im.get_cmap()
            current_cmap.set_bad(color = 'white')
            cbar.ax.tick_params(labelsize = 12)
            ylabels = [str(i) + '\%' for i in list(np.arange(-4, 5)*2.5)]
            cbar.ax.set_yticklabels(ylabels)
        else:
            cbar.remove()
        # Add percent deviation inside each box
        for i in range(len(discrp)):
            for j in range(len(discrp)):
                if not np.isnan(discrp[i, j]):
                    text = str('{:.1f}'.format(round(discrp[i, j], 2)))
                    ax.text(j, i, text, color = 'black',
                            ha = 'center', va = 'center')
        # Add formula
        text = (r'$\frac{\mathrm{' + a + '}-\mathrm{' + b + '}}'
                + '{\mathrm{' + b + '}} (\%)$')
        ax.text(4.7, 0.1, text, fontsize = 12, color = 'black',
                ha = 'center', va = 'center')
        # Add Battleship-style coordinates as ticks
        ax.xaxis.tick_bottom()
        xlabels = ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']
        ylabels = list(np.arange(1, 9))
        plt.setp(ax, xticks = np.arange(8), xticklabels = xlabels,
                 yticks = np.arange(8), yticklabels = ylabels)
        ax.tick_params(axis = 'both', which = 'both', length = 0)
        #---
        #  Add a title
        #---
        if lang == 'en':
            fig.suptitle('Relative difference in power distribution between\n'
                         + a + ' and ' + b + ' on Tihange\nfirst zero-power '
                         'start-up, ' + CRtext[controlrod])
        #---
        #  Save plot as pdf (vectorized)
        #---
        fig.savefig('output_PowDifferenceTih/' + controlrod + '_' + a + '_vs_'
                    + b + '.pdf', bbox_inches='tight')
        #---
        #  Clean-up for next plot
        #---
        plt.close('all')
        del fig
        del ax
        # Compute and print the maximum relative deviation between the two
        # selected sources
        diffmax = np.max(np.abs(discrp1D))
        print(a + '/' + b, '{:02.2f}'.format(diffmax))
