#
#  Plotting power relative difference on 2D maps
#  Usage : python3 PowDifference.py
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
import itertools
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from BasicFunctions import *
import serpentTools
from serpentTools.settings import rc
#---
#  Initialize the container of all combinations of power distributions
#---
controlrods = ['ARO', 'D', 'CD']
CRtext = {}
CRtext['ARO'] = 'all rods out'
CRtext['D'] = 'D rod bank inserted'
CRtext['CD'] = 'C and D rod banks inserted'
sources = ['Drakkar', 'Serpent', 'Framatome', 'EDF']
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
for file in list(glob.glob('../Serpent/FullCore/*.sss2_det0.m')):
    cbor = file.split('/')[-1].split('_')[0]
    cbor = int(re.findall(r'\d+', cbor)[0])
    if cbor == 1206:
        controlrod = 'ARO'
    elif cbor == 1084:
        controlrod = 'D'
    elif cbor == 960:
        controlrod = 'CD'
    # Retrieve Serpent detector output
    det = serpentTools.read(file)
    # Total energy deposition, see
    # http://serpent.vtt.fi/mediawiki/index.php/ENDF_reaction_MT%27s_and_macroscopic_reaction_numbers
    reaction = '-80'
    power = det.detectors[reaction].tallies
    if reaction == '-80':
        power = power[0]
    # Number of assembly-sized tallies outside the core, on each side
    outside = 2
    # Produce the core layout corresponding to tallies
    hstack = np.zeros((len(FullCoreLayout), outside))
    TalliesLayout = np.hstack([hstack, FullCoreLayout, hstack])
    vstack = np.zeros((outside, len(FullCoreLayout) + 2*outside))
    TalliesLayout = np.vstack([vstack, TalliesLayout, vstack])
    # Compute fraction of power delivered in the reflector
    corepower = np.where(TalliesLayout == 0, 0, power[1])
    radialreflpower = np.where(TalliesLayout == 1, 0, power[1])
    botreflpower = np.sum(power[0])
    topreflpower = np.sum(power[2])
    reflpower = np.sum(radialreflpower) + botreflpower + topreflpower
    print('Fraction of the power Serpent delivers in the reflector = '
          + str(np.sum(reflpower)/np.sum(power)) + ' (' + controlrod + ')')
    # Ignore reflector power and renormalize core power
    corepower = corepower/np.sum(corepower)*np.count_nonzero(corepower)
    # The input is symmetrical by eighth, so we can symmetrize (a fold followed
    # by an unfold) in order to increase the statistical strength
    symcorepower = UnfoldEighth(FoldEighth(corepower))
    # Estimate the uncertainty
    sym = symcorepower
    unsym = corepower
    reldif = unsym[np.nonzero(unsym)] / sym[np.nonzero(sym)] - 1
    print("min, max (%)=" + str(np.min(reldif)*100) + ", "
          + str(np.max(reldif)*100))
    # Remove all powers equal to zero
    powers[controlrod]['Serpent'] = symcorepower[np.nonzero(symcorepower)]
#---
#  Drakkar power distributions
#---
for controlrod in controlrods:
    os.system('ln -s ../Drakkar/Linux_x86_64/_Power' + controlrod
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
#  Plot all the power distributions
#---
os.system('mkdir -p output_PowDifference')
pmin = math.inf
pmax = -math.inf
for controlrod in controlrods:
    for source in sources:
        pmin = min(pmin, np.min(powers[controlrod][source]))
        pmax = max(pmax, np.max(powers[controlrod][source]))
for controlrod in controlrods:
    for source in sources:
        # Show only a quarter of the core
        powerplot = fold(Deploy2D(powers[controlrod][source], FullCoreLayout))
        # Replace every zeroes with NaNs, in order to easily get them colored
        # in white (with set_bad)
        powerplot = np.where(np.isclose(powerplot, 0), float('nan'), powerplot)
        # Produce the plot
        fig, ax = plt.subplots(1, figsize=[6, 6])
        axins = inset_axes(ax,
                           width = "3%",
                           height = "100%",
                           loc = 3,
                           bbox_to_anchor = (1.02, 0., 1, 1),
                           bbox_transform = ax.transAxes,
                           borderpad = 0,
                           )
        im = ax.imshow(powerplot, interpolation = 'nearest', cmap = 'coolwarm',
                       vmin = pmin, vmax = pmax)
        cbar = plt.colorbar(im, cax = axins, orientation = 'vertical')
        current_cmap = im.get_cmap()
        current_cmap.set_bad(color = 'white')
        cbar.ax.tick_params(labelsize = 15)
        # Add power in each assembly
        fontsize_ = 10
        for i in range(len(powerplot)):
            for j in range(len(powerplot)):
                if not np.isnan(powerplot[i, j]):
                    text = str('{:.2f}'.format(round(powerplot[i, j], 2)))
                    ax.text(j, i, text, fontsize = fontsize_, color = 'black',
                            ha = 'center', va = 'center')
        # Add Battleship-style coordinates as ticks
        ax.xaxis.tick_bottom()
        xlabels = ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']
        ylabels = list(np.arange(1, 9))
        plt.setp(ax, xticks = np.arange(8), xticklabels = xlabels,
                 yticks = np.arange(8), yticklabels = ylabels)
        #---
        #  Add a title
        #---
        if source == 'EDF' or source == 'Framatome':
            txt = source + ' pseudo-measured'
        else:
            txt = source + ' JEFF-3.3'
        fig.suptitle(txt + ' power distribution on Tihange\nFirst zero-power '
                     + 'start-up, ' + CRtext[controlrod])
        #---
        #  Save plot as pdf (vectorized)
        #---
        fig.savefig('output_PowDifference/' + source + '_' + controlrod
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
        fig, ax = plt.subplots(1, figsize=[6, 6])
        axins = inset_axes(ax,
                           width = "3%",
                           height = "100%",
                           loc = 3,
                           bbox_to_anchor = (1.02, 0., 1, 1),
                           bbox_transform = ax.transAxes,
                           borderpad = 0,
                           )
        im = ax.imshow(discrp, interpolation = 'nearest', cmap = 'coolwarm',
                       vmin = -15, vmax = 15)
        cbar = plt.colorbar(im, cax = axins, orientation = 'vertical')
        current_cmap = im.get_cmap()
        current_cmap.set_bad(color = 'white')
        cbar.ax.tick_params(labelsize = 15)
        # Add percent deviation inside each box
        fontsize_ = 10
        for i in range(len(discrp)):
            for j in range(len(discrp)):
                if not np.isnan(discrp[i, j]):
                    text = str('{:.2f}'.format(round(discrp[i, j], 2)))
                    ax.text(j, i, text, fontsize = fontsize_, color = 'black',
                            ha = 'center', va = 'center')
        # Add formula
        text = (r'$\frac{\mathrm{' + a + '}-\mathrm{' + b + '}}'
                + '{\mathrm{' + b + '}} (\%)$')
        ax.text(5.5, 0, text, fontsize = fontsize_ + 6, color = 'black',
                ha = 'center', va = 'center')
        # Add Battleship-style coordinates as ticks
        ax.xaxis.tick_bottom()
        xlabels = ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']
        ylabels = list(np.arange(1, 9))
        plt.setp(ax, xticks = np.arange(8), xticklabels = xlabels,
                 yticks = np.arange(8), yticklabels = ylabels)
        #---
        #  Add a title
        #---
        fig.suptitle('Relative difference in power distribution between\n'
                     + a + ' and ' + b + ' on Tihange\nFirst zero-power '
                     'start-up, ' + CRtext[controlrod])
        #---
        #  Save plot as pdf (vectorized)
        #---
        fig.savefig('output_PowDifference/' + controlrod + '_' + a + '_vs_' + b
                    + '.pdf', bbox_inches='tight')
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
