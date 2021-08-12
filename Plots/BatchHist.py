#
#  Plot batch history, to visualize its convergence
#  Usage : python3 BatchHist.py
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import numpy as np
import os
import math
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from BasicFunctions import *
import re
import serpentTools
from serpentTools.settings import rc
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en
#
nbatchplot = 198
#
controlrods = ['CD', 'D', 'ARO']
text = {}
text['ARO'] = 'all rods out'
text['D'] = 'D rod bank inserted'
text['CD'] = 'C and D rod banks inserted'
# Storage for independant runs' results
indruns = {}
for controlrod in controlrods:
    indruns[controlrod] = []
# Layout of the assemblies in the core
CoreLayout, FullCoreLayout = GetCoreLayout(157)
SerpentLayout = GetSerpentLayout(FullCoreLayout, 2)
for file in glob.glob('../Serpent/BatchHist/Tihange*ppm_*.sss2_his0.m'):
    print(file)
    # Identify the file we're about to read
    basename = file.split('/')[-1].split('.')[0]
    cbor = int(re.sub("[^0-9]", "", basename.split('_')[0]))
    if cbor == 1206:
        controlrod = 'ARO'
    elif cbor == 1084:
        controlrod = 'D'
    elif cbor == 960:
        controlrod = 'CD'
    else:
        raise Exception(cbor + ' unknown.')
    # Read that Serpent history file
    hist = serpentTools.read(file)
    # For a detector (here -80), the results are grouped in three columns that
    # provide (1) the cycle-wise value,
    #         (2) the cumulative mean and
    #         (3) the corresponding relative statistical error.
    # Cf. http://serpent.vtt.fi/mediawiki/index.php/Description_of_output_files#History_output
    # Only the first is kept here.
    cycle = hist['det-80'][:, 0::3]
    # Among that, only the total deposited energy is kept
    cycle = cycle[:, 0::3]
    # Among that, only the central layer (the active core) is kept ; top and
    # bottom layers are splitted out
    cycle = np.split(cycle, 3, axis = 1)[1]
    # Raise exception if inactive generations are being used. One (1) inactive
    # generation is the minimum permitted by Serpent and therefore the only
    # one we accept.
    if np.isclose(0.0, np.sum(cycle[1]), atol = 1e-20):
        raise Exception('The convergence plot cannot be achieved since '
                        + 'inactive generations have been detected: '
                        + '*set pop [x] [y] 1* should be used in Serpent.')
    # Reshaping into nbatch*size*size
    nbatch = np.shape(cycle)[0]
    nxysize = int(math.sqrt(np.shape(cycle)[1]))
    cycle.shape = (nbatch, nxysize, nxysize)
    # Prepare quarter of cores
    nxyassembly = np.shape(FullCoreLayout)[0]
    cycle4th = np.zeros((nbatch, int((nxyassembly - 1)/2) + 1,
                         int((nxyassembly - 1)/2) + 1))
    for ibatch in range(1, nbatch):
        # Filter out the reflector to get only the response of the core
        thiscycle = np.where(SerpentLayout == 0, 0, cycle[ibatch])
        # Normalize the core power
        nassembly = np.count_nonzero(FullCoreLayout)
        thiscycle = nassembly * thiscycle / np.sum(thiscycle)
        # Remove columns and rows full of zeroes
        thiscycle = Deploy2D(thiscycle[np.nonzero(thiscycle)], FullCoreLayout)
        # The input is symmetrical by eighth, so we can symmetrize (a fold
        # followed by an unfold) in order to increase the statistical strength
        thiscycle = UnfoldEighth(FoldEighth(thiscycle))
        # For plots, keep only a quarter of the result
        cycle4th[ibatch] = fold(thiscycle)
    # Store the results of this file
    indruns[controlrod].append(cycle4th)
tmp = []
for controlrod in controlrods:
    # Convert this list of n-D NumPy array into a (n+1)-D NumPy array
    indruns[controlrod] = np.array(indruns[controlrod])
    # Update 'controlrods' list so that only those states that were indeed used
    # can be the subject of loops to come
    if np.size(indruns[controlrod]) != 0:
        tmp.append(controlrod)
controlrods = tmp
for controlrod in controlrods:
    # Mean over independant runs
    indmean = np.mean(indruns[controlrod], axis = 0)
    # Mean over active generations
    indmeanact = np.mean(indmean[nbatchplot:], axis = 0)
    #---
    #  Plotting
    #---
    # Initialize figure
    fig, axs = plt.subplots(8, 8, sharex = 'all',
                            figsize = set_size('square'),
                            gridspec_kw = {'hspace': 0, 'wspace': 0})
    xaxis = np.arange(nbatchplot - 1)
    for x in range(0, len(CoreLayout[0, :])):
        for y in range(0, len(CoreLayout[:, 0])):
            if CoreLayout[x, y] == 0 and [x, y] != [1, 7]:
                axs[x, y].set_axis_off()
            else:
                # Add a black vertical line on zero
                axs[x, y].axhline(y = 0, linestyle = '--', linewidth = 0.5,
                                  color = 'k')
                if [x, y] != [1, 7]:
                    # Compute relative discrepancy with regard to mean over
                    # active generations
                    reldiscr = ((indmean[1:nbatchplot, x, y]/indmeanact[x, y]
                                 - 1)*100)
                else:
                    # Add an outside-of-core example for axis labels
                    reldiscr = ((indmean[1:nbatchplot, 2, 2]/indmeanact[2, 2]
                                 - 1)*100)
                axs[x, y].plot(xaxis, reldiscr)
                # Rotate labels
                axs[x, y].tick_params(axis = 'x', labelrotation = 90)
    for x in range(0, len(CoreLayout[0, :])):
        for y in range(0, len(CoreLayout[:, 0])):
            # Adjust x-axis ticks
            axs[x, y].set_xlim(0, nbatchplot + 20)
            plt.setp(axs, xticks=[0, 75, 150])
            # Adjust y-axis ticks
            axs[x, y].set_ylim(-6.5, +6.5)
            if y == 0:
                labels = axs[x, y].get_yticks()
                # Add percent sign after yticks
                labels = [str(i) + '\%' for i in labels]
                # Prefix with a plus sign for labels not starting with a minus
                # sign
                labels = ['+' + i if i[0] != '-' else i for i in labels]
                # Add dollar signs, so that math TeX typography is used for
                # plus and minus signs
                labels = ['$' + i + '$' for i in labels]
                # For central value, keep it simple : 0 (no sign, no percent)
                labels[int(len(labels)/2 - 0.5)] = '0'
                axs[x, y].set_yticklabels(labels)
            else:
                axs[x, y].set_yticks([])
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
    #  Add the axis labels on an outside-of-core example
    #---
    axs[1, 7].tick_params(axis = 'x', labelbottom = True)
    if lang == 'en':
        labelx = 'Inactive\ngeneration\nindex'
        labely = ('Relative error on assembly\n'
                  + 'power with respect to\n'
                  + '1000 active generations')
    elif lang == 'fr':
        labelx = 'Indice de\ngénération\ninactive'
        labely = ('Erreur relative sur la\n'
                  + 'puissance de l\'assemblage,\n'
                  + 'par rapport aux 1000\n'
                  + 'générations actives')
    axs[1, 7].set_xlabel(labelx)
    axs[1, 7].set_ylabel(labely)
    #---
    #  Add a title and save plot as pdf (vectorized)
    #---
    if lang == 'en':
        title = ('Convergence of assembly power during 200 inactive '
                 + 'generations,'
                 + '\nuniformly initialized, ' + r'$10^5$ neutrons/generation,'
                 + ' averaged over\n'
                 + '50 independant simulations with different random seeds\n'
                 + 'on Tihange first zero-power start-up, ' + text[controlrod])
    elif lang == 'fr':
        title = ''
    os.system('mkdir -p output_BatchHist')
    if title:
        st = fig.suptitle(title, y = 1.01)
        fig.savefig('output_BatchHist/' + controlrod + '.pdf',
                    extra_artists = [st], bbox_inches = 'tight')
    else:
        fig.savefig('output_BatchHist/' + controlrod + '.pdf',
                    bbox_inches = 'tight')
    #---
    #  Clean-up for next plot
    #---
    plt.close('all')
print("Plotting completed")
