#
#  Plotting power on 2D maps for TMC samples (both Drakkar and Serpent)
#  Usage : python3 TMC_comparison.py
#  Author ; V. Salino (IRSN), 06/2021
#

# Imports
import lcm
import numpy as np
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import ks_2samp
import math
import os
import glob
import re
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from BasicFunctions import *
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en

# Initialize a gaussian vector (used for legends)
gaussian = np.random.randn(300)
# Declare control rod insertions and its full names
controlrods = ['CD', 'D', 'ARO']
text = {}
text['ARO'] = 'all rods out'
text['D'] = 'D rod bank inserted'
text['CD'] = 'C and D rod banks inserted'
for controlrod in controlrods:
    print('controlrod = ' + controlrod)
    if controlrod == 'ARO':
        cbor = 1206
    elif controlrod == 'D':
        cbor = 1084
    elif controlrod == 'CD':
        cbor = 960
    #---
    #  Drakkar power distributions
    #---
    # Create links toward every available Drakkar power distribution
    os.system('ln -s ../Drakkar/Output_TMC/_Power' + controlrod + '_*.ascii .')
    # List isotopes we've been (potentially) randomly sampling in Drakkar
    firstfile = glob.glob('_Power' + controlrod + '_*.ascii')[0]
    ResultFile = lcm.new('LCM_INP', firstfile[1:])
    drkisotopes = ResultFile['NamIso'].split()
    del ResultFile
    # Initialize dictionnary (isotope-wise) of lists (sampling-wise) that will
    # contain the retrieved Drakkar data
    powers = {}
    for iso in drkisotopes:
        powers[iso] = []
    for file in list(glob.glob('_Power' + controlrod + '_*.ascii')):
        # Loading data while removing beginning-of-filename underscore
        ResultFile = lcm.new('LCM_INP', file[1:])
        os.system("rm " + file)
        # Compute mean in order to normalize power to 1.0 per average assembly
        power = ResultFile['POWER-CHAN']
        mean = np.sum(power)/len(power)
        power = power/mean
        # Determine which isotope was randomly sampled in order to store the
        # obtained power distribution in the right place
        for iso in drkisotopes:
            if (int(ResultFile[iso]) != -33):
                powers[iso].append(power)
        del ResultFile
    # Rename 'H1_H2O' isotope into 'H1', so that its name in Drakkar and
    # Serpent matches
    if 'H1_H2O' in drkisotopes:
        drkisotopes.remove('H1_H2O')
        drkisotopes.append('H1')
        powers['H1'] = powers.pop('H1_H2O')
    # Convert this list of 1D numpy array into a 2D numpy array:
    # [irand, iassembly]
    tmp = {}
    for iso in drkisotopes:
        tmp[iso] = np.array(powers[iso])
    powers = tmp
    # Update 'isotopes' list so that only those isotopes that were indeed
    # randomly sampled can be the subject of loops to come
    isotopes_used = []
    for iso in drkisotopes:
        if powers[iso].size != 0:
            isotopes_used.append(iso)
    drkisotopes = isotopes_used
    # Layout of the assemblies in the core
    nbassembly = len(powers[drkisotopes[0]][0, :])
    CoreLayout, FullCoreLayout = GetCoreLayout(nbassembly)
    #
    FullPowers2D = {}
    Powers2D = {}
    for iso in drkisotopes:
        # Deploy the 1D-indexed power maps into 2D-indexed power maps
        FullPowers2D[iso] = np.zeros(( len(powers[iso][:, 0]),      # irand
                                       len(FullCoreLayout[:, 0]),   # Y-axis
                                       len(FullCoreLayout[0, :]) )) # X-axis
        i = 0
        for x in range(0, len(FullCoreLayout[0, :])):
            for y in range(0, len(FullCoreLayout[:, 0])):
                if FullCoreLayout[x, y] != 0:
                    FullPowers2D[iso][:, x, y] = powers[iso][:, i]
                    i = i + 1
        # Fold in quarter maps (averaging per rotation) for plotting purposes
        Powers2D[iso] = np.zeros(( len(powers[iso][:, 0]),  # irand
                                   len(CoreLayout[:, 0]),   # Y-axis
                                   len(CoreLayout[0, :]) )) # X-axis
        for irand in range(0, len(FullPowers2D[iso][:, 0, 0])):
            Powers2D[iso][irand, :, :] = fold(FullPowers2D[iso][irand, :, :])
    #---
    #  Serpent power distributions
    #---
    # Initialize dictionnary (isotope-wise) of lists (sampling-wise) that will
    # contain the retrieved data
    powers_sss = {}
    # List the isotopes that have been the subject of a TMC within Serpent
    sssisotopes = []
    for file in list(glob.glob('../Serpent/TMC/TihangeFullCore' + str(cbor)
                               + 'ppm_*.sss2_det0.m')):
        iso = file.split('/')[-1].split('_')[2]
        if iso not in sssisotopes:
            sssisotopes.append(iso)
            powers_sss[iso] = []
        power = ReadPdistrSerpent(file, FullCoreLayout, controlrod)
        powers_sss[iso].append(fold(Deploy2D(power, FullCoreLayout)))
    # Convert this list of 1D numpy array into a 2D numpy array:
    # [irand, iassembly]
    tmp = {}
    for iso in sssisotopes:
        tmp[iso] = np.array(powers_sss[iso])
    powers_sss = tmp
    # Check that each isotope that has been the subject of a TMC within
    # Serpent is also available with Drakkar
    for sssiso in sssisotopes:
        if sssiso not in drkisotopes:
            raise Exception (sssiso + 'has been sampled in Serpent, but does '
                             + 'not seem to be available for Drakkar. The '
                             + 'isotopes available for Drakkar are:'
                             + drkisotopes)
    #---
    #  Plotting
    #---
    for iso in sssisotopes:
        print(iso)
        # Initialize figure
        fig, axs = plt.subplots(8, 8, sharex = 'all',
                                figsize = set_size('square'),
                                gridspec_kw = {'hspace': 0, 'wspace': 0})
        #---
        #  Add the data
        #---
        for x in range(0, len(CoreLayout[0, :])):
            for y in range(0, len(CoreLayout[:, 0])):
                if CoreLayout[x, y] == 0:
                    axs[x, y].set_axis_off()
                else:
                    # Plot an histogram with the Drakkar assembly power
                    axs[x, y].hist(Powers2D[iso][:, x, y], bins = 20,
                                   density = True)
                    # Compute the Drakkar relative standard deviation for each
                    # assembly
                    mean = np.mean(Powers2D[iso][:, x, y])
                    relstd = np.std(Powers2D[iso][:, x, y], ddof = 1)/mean*100
                    textstr = r'$\sigma_D$=$%.1f$' % (relstd, ) + '\%'
                    # Plot an histogram with the Serpent assembly power
                    axs[x, y].hist(powers_sss[iso][:, x, y], bins = 20,
                                   density = True)
                    # Compute the Serpent relative standard deviation for each
                    # assembly
                    mean = np.mean(powers_sss[iso][:, x, y])
                    relstd = np.std(powers_sss[iso][:, x, y], ddof = 1)/mean*100
                    textstr += '\n' + r'$\sigma_S$=$%.1f$' % (relstd, ) + '\%'
                    # Compute the Kolmogorov-Smirnov statistic and its p-value
                    # on Serpent and Drakkar samples. We know that Serpent mean
                    # is different than Drakkar mean. We rather want to compare
                    # the higher orders. The means of each one are thus
                    # extracted, before computing the Kolmogorov-Smirnov
                    # statistic.
                    D, p = ks_2samp(powers_sss[iso][:, x, y] -
                                     np.mean(powers_sss[iso][:, x, y]),
                                     Powers2D[iso][:, x, y] -
                                     np.mean(Powers2D[iso][:, x, y]))
                    textstr += '\n$D$=$%.2f$' % (D)
                    textstr += '\n$p$=$%.2f$' % (p)
                    # Print statistics for each assembly
                    axs[x, y].text(0.05, 0.92, textstr, fontsize = 10,
                                   transform=axs[x, y].transAxes,
                                   verticalalignment='top')
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
        #---
        #  Add a legend
        #---
        if lang == 'en':
            textstr = (r'$\sigma_D$ : Drakkar relative standard deviation'
                       + ' (in \%)\n'
                       + r'$\sigma_S$ : Serpent relative standard deviation'
                       + ' (in \%)\n'
                       + '$D$ : Kolmogorov-Smirnov statistic' + '\n'
                       + '$p$ : Kolmogorov-Smirnov p-value')
            axs[0, 4].text(-0.3, 0.95, textstr,
                           transform = axs[0, 4].transAxes,
                           verticalalignment = 'top', fontsize = 10)
        elif lang == 'fr':
            textstr = (r'$\sigma_D$ : écart-type relatif de Drakkar'
                       + ' (en \%)\n'
                       + r'$\sigma_S$ : écart-type relatif de Serpent'
                       + ' (en \%)\n'
                       + '$D$ : statistique de Kolmogorov-Smirnov' + '\n'
                       + '$p$ : valeur-p de Kolmogorov-Smirnov')
            axs[0, 2].text(0.05, 1.04, textstr,
                           transform = axs[0, 2].transAxes,
                           verticalalignment = 'top', fontsize = 10)
        #---
        #  Add the axis labels on an fictious (gaussian) distribution example
        #---
        axs[1, 6].set_axis_on()
        axs[1, 6].hist(1.0225 + 0.015 * gaussian, bins = 20, density = True)
        axs[1, 6].set_yticks([])
        axs[1, 6].set_aspect(1./axs[1, 6].get_data_ratio())
        if lang == 'en':
            axs[1, 6].set_xlabel('Normalized assembly\npower')
            axs[1, 6].set_ylabel('Probability\ndensity')
        elif lang == 'fr':
            axs[1, 6].set_xlabel('Puissance normalisée\nde l\'assemblage')
            axs[1, 6].set_ylabel('Densité\nde\nprobabilité')
        axs[1, 6].tick_params(axis = 'x', labelbottom = True)
        #---
        #  Add a legend
        #---
        label1 = matplotlib.lines.Line2D([], [], color = '#1f77b4',
                                         label = 'Drakkar')
        label2 = matplotlib.lines.Line2D([], [], color = '#ff7f0e',
                                         label = 'Serpent')
        axs[0, 7].legend(handles=[label1, label2])
        #---
        #  Add a title and save plot as pdf (vectorized)
        #---
        os.system('mkdir -p output_TMC_comparison')
        if lang == 'en':
            title = ('Probability density vs. normalized assembly power, '
                     + 'with Drakkar, Serpent\nand '
                     + 'JEFF-3.3 best estimate except for '
                     + iso.replace('_', '\_') + ' ('
                     + str(len(Powers2D[iso][:, 0, 0]))
                     + ' random samples)\non Tihange first zero-power '
                     + 'start-up, ' + text[controlrod])
            st = fig.suptitle(title)
            fig.savefig('output_TMC_comparison/' + controlrod + '_' + iso
                        + '.pdf', bbox_inches='tight', bbox_extra_artists=[st])
        elif lang == 'fr':
            fig.savefig('output_TMC_comparison/' + controlrod + '_' + iso
                        + '.pdf', bbox_inches='tight')
        #---
        #  Clean-up for next plot
        #---
        plt.close('all')
    #---
    #  Merge pdf with ghostscript
    #---
    os.system('gs -q -dNOPAUSE -sDEVICE=pdfwrite '
              + '-sOUTPUTFILE=output_TMC_comparison/' + controlrod + '.pdf '
              + '-dBATCH output_TMC_comparison/' + controlrod + '_*.pdf')
