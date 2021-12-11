#
#  Plotting power relative difference on 2D maps, for CP0 cases
#  Usage : python3 PowDifferenceCP0.py
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

# Loop over the two cases we're interested in
cases = ['Almaraz', 'Fessenheim-Bugey']
sources = {}
powers = {}
for case in cases:
    #---
    #  Initialize the container of all combinations of power distributions
    #---
    if case == 'Almaraz':
        sources[case] = ['Serpent', 'Mesures']
    elif case == 'Fessenheim-Bugey':
        sources[case] = ['Drakkar', 'Serpent']
    powers[case] = {}
    if case == 'Almaraz':
        #---
        #  Measurement-derived power distributions
        #---
        file = '../Measurements/CP0/Almaraz-2_Cycle1_0MWdt.csv'
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
        powers[case]['Mesures'] = meas
        CoreLayout, FullCoreLayout = GetCoreLayout(len(meas[meas>0]))
        FullCoreLayout4X = np.repeat(np.repeat(FullCoreLayout,
                                               2, axis = 0), 2, axis = 1)
    #---
    #  Serpent power distributions
    #---
    links = glob.glob('../Serpent/SuperBatch/' + case
                      + 'FullCore*ppm_*.sss2_det0.m')
    # Storage for independant runs' results
    indruns = []
    for link in links:
        indruns.append(ReadPdistrSerpent(link, FullCoreLayout))
    powers[case]['Serpent'] = Deploy2D(np.mean(indruns, axis = 0),
                                       FullCoreLayout)
    if case == 'Fessenheim-Bugey':
        #---
        #  Drakkar power distributions
        #---
        os.system('ln -s ../Drakkar/Output_FSH-BUG_BestEstimate/_PowerARO'
                  + '_*.ascii .')
        for file in list(glob.glob('_PowerARO_*.ascii')):
            # Loading data while removing beginning-of-filename underscore
            ResultFile = lcm.new('LCM_INP', file[1:])
            os.system("rm " + file)
            power = ResultFile['POWER-CHAN']
            # Compute mean in order to normalize power to 1.0 per average assembly
            mean = np.sum(power)/len(power)
            power = power/mean
            # Each FullCoreLayout value is splitted in 4 values
            power = matrixaverage(Deploy2D(power, FullCoreLayout4X), 2)/2
            powers[case]['Drakkar'] = power
#---
#  Plot all the power distributions
#---
os.system('mkdir -p output_PowDifferenceCP0')
pmin = math.inf
pmax = -math.inf
for case in cases:
    for source in sources[case]:
        power = powers[case][source][powers[case][source]>0]
        pmin = min(pmin, np.min(power))
        pmax = max(pmax, np.max(power))
for case in cases:
    for source in sources[case]:
        plot = fold(powers[case][source])
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
                           borderpad = 0,
                           )
        im = ax.imshow(plot, interpolation = 'nearest', cmap = 'coolwarm',
                       vmin = pmin, vmax = pmax)
        cbar = plt.colorbar(im, cax = axins, orientation = 'vertical')
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
            if source == 'Mesures':
                txt = 'Pseudo-measured'
            else:
                txt = source + ' JEFF-3.3'
            txt = (txt + ' power distribution on ' + case
                   + '\nfirst zero-power start-up')
            if source == 'Serpent':
                nbruns = np.shape(indruns)[0]
                txt += ('\n' + r'Relative standard error of the mean '
                        + '($1\sigma$) from ' + str(nbruns) + '\n'
                        + 'independant simulations with different random seeds')
                plt.subplots_adjust(top = 0.85)
            st = fig.suptitle(txt)
            fig.savefig('output_PowDifferenceCP0/' + case + '_' + source
                        + '.pdf', bbox_inches='tight', extra_artists=[st])
        elif lang == 'fr':
            fig.savefig('output_PowDifferenceCP0/' + case + '_' + source
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
    # Get every two-by-two combinations between every sources
    for (a, b) in list(itertools.combinations(sources[case], 2)):
        # Compute relative deviation (discrepancy) between the two selected
        # sources ('a' and 'b')
        discrp = (powers[case][a]/powers[case][b] - 1)*100
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
                           borderpad = 0,
                           )
        im = ax.imshow(discrp, interpolation = 'nearest', cmap = 'coolwarm',
                       vmin = -10, vmax = 10)
        cbar = plt.colorbar(im, cax = axins, orientation = 'vertical')
        current_cmap = im.get_cmap()
        current_cmap.set_bad(color = 'white')
        cbar.ax.tick_params(labelsize = 12)
        ylabels = [str(i) + '\%' for i in list(np.arange(-4, 5)*2.5)]
        cbar.ax.set_yticklabels(ylabels)
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
                         + a + ' and ' + b + ' on ' + case
                         + '\nfirst zero-power start-up')
        #---
        #  Save plot as pdf (vectorized)
        #---
        fig.savefig('output_PowDifferenceCP0/' + case + '_' + a + '_vs_' + b
                    + '.pdf', bbox_inches='tight')
        #---
        #  Clean-up for next plot
        #---
        plt.close('all')
        del fig
        del ax
        # Compute and print the maximum relative deviation between the two
        # selected sources
        diffmax = np.max(np.abs(discrp[~np.isnan(discrp)]))
        print(a + '/' + b, '{:02.2f}'.format(diffmax))
