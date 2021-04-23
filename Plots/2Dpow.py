#
#  Plotting power uncertainty on 2D maps
#  Usage : python3 2Dpow.py
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import lcm
import numpy as np
from scipy.stats import skew
from scipy.stats import kurtosis
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from BasicFunctions import *
# Initialize a gaussian vector (used for legends)
gaussian = np.random.randn(300)
# Declare control rod insertions and its full names
controlrods = ['CD', 'D', 'ARO']
text = {}
text['ARO'] = 'all rods out'
text['D'] = 'D rod bank inserted'
text['CD'] = 'C and D rod banks inserted'
# Create a dictionary to store the standard deviation of the assembly with the
# most uncertain power, among all those assemblies in the core.
maxstd = {}
for controlrod in controlrods:
    print('controlrod = ' + controlrod)
    maxstd[controlrod] = {}
    # Create links toward every available power distribution
    os.system('ln -s ../Drakkar/Linux_x86_64/_Power' + controlrod
              + '_*.ascii .')
    # List isotopes we've been (potentially) randomly sampling
    firstfile = glob.glob('_Power' + controlrod + '_*.ascii')[0]
    ResultFile = lcm.new('LCM_INP', firstfile[1:])
    isotopes = ResultFile['NamIso'].split()
    del ResultFile
    # Initialize dictionnary (isotope-wise) of lists (sampling-wise) that will
    # contain the retrieved data
    powers = {}
    for iso in isotopes:
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
        for iso in isotopes:
            if int(ResultFile[iso]) != -311:
                powers[iso].append(power)
        del ResultFile
    # Convert this list of 1D numpy array into a 2D numpy array:
    # [irand, iassembly]
    tmp = {}
    for iso in isotopes:
        tmp[iso] = np.array(powers[iso])
    powers = tmp
    # Update 'isotopes' list so that only those isotopes that were indeed
    # randomly sampled can be the subject of loops to come
    isotopes_used = []
    for iso in isotopes:
        if powers[iso].size != 0:
            isotopes_used.append(iso)
    isotopes = isotopes_used
    # Layout of the assemblies in the core
    CoreLayout, FullCoreLayout = GetCoreLayout(len(powers[isotopes[0]][0, :]))
    #
    FullPowers2D = {}
    Powers2D = {}
    for iso in isotopes:
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
    #  Plotting
    #---
    for iso in isotopes:
        print(iso)
        # Initialize maximum relative standard deviation to zero
        maxstd[controlrod][iso] = 0
        # Initialize figure
        fig, axs = plt.subplots(8, 8, sharex = 'all', figsize = (8, 8),
                                gridspec_kw = {'hspace': 0, 'wspace': 0})
        #---
        #  Add the data
        #---
        for x in range(0, len(CoreLayout[0, :])):
            for y in range(0, len(CoreLayout[:, 0])):
                if CoreLayout[x, y] == 0:
                    axs[x, y].set_axis_off()
                else:
                    # Plot an histogram with the assembly power
                    axs[x, y].hist(Powers2D[iso][:, x, y], bins = 20,
                                   density = True)
                    # Compute relative standard deviation for each assembly
                    mean = np.mean(Powers2D[iso][:, x, y])
                    relstd = np.std(Powers2D[iso][:, x, y])/mean*100
                    textstr = (r'$\mu$=%.2f' % mean + '\n'
                               + r'$\sigma$=%.1f%%' % (relstd, ))
                    # Store its maximum for future purposes (summary plot)
                    maxstd[controlrod][iso] = max(relstd,
                                                  maxstd[controlrod][iso])
                    # Compute skewness and kurtosis only if sigma is large
                    # enough (> 0.1%)
                    if relstd > 0.1:
                        skwn = skew(Powers2D[iso][:, x, y])
                        kurt = kurtosis(Powers2D[iso][:, x, y])
                        textstr +=  '\nS=%.1f' % skwn + '\nK=%.1f' % kurt
                    # Print statistics for each assembly
                    axs[x, y].text(0.05, 0.95, textstr,
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
        textstr = (r'$\mu$ : Mean of the power assembly' + '\n'
                   + r'$\sigma$ : Relative standard deviation (in %)' + '\n'
                   + r'S : Skewness (shown if $\sigma > 0.1$%)' + '\n'
                   + r'K : Excess kurtosis (shown if $\sigma > 0.1$%)')
        axs[0, 4].text(0.6, 0.95, textstr,
                       transform = axs[0, 4].transAxes,
                       verticalalignment = 'top')
        #---
        #  Add the axis labels on an fictious (gaussian) distribution example
        #---
        axs[2, 7].set_axis_on()
        axs[2, 7].hist(1.0225 + 0.015 * gaussian, bins = 20, density = True)
        axs[2, 7].set_yticks([])
        axs[2, 7].set_aspect(1./axs[2, 7].get_data_ratio())
        axs[2, 7].set_xlabel('Normalized\nassembly\npower')
        axs[2, 7].set_ylabel('Probability\ndensity')
        axs[2, 7].tick_params(axis = 'x', labelbottom = True)
        #---
        #  Add a title
        #---
        st = fig.suptitle('Probability density vs. normalized assembly power, '
                          + 'with Drakkar and \n'
                          + 'JEFF-3.3 best estimate except for ' + iso + ' ('
                          + str(len(Powers2D[iso][:, 0, 0]))
                          + ' random samples) on\nTihange first zero-power '
                          + 'start-up, ' + text[controlrod])
        #---
        #  Save plot as pdf (vectorized)
        #---
        os.system('mkdir -p output_2Dpow')
        fig.savefig('output_2Dpow/' + controlrod + '_' + iso + '.pdf',
                    bbox_extra_artists=[st], bbox_inches='tight')
        #---
        #  Clean-up for next plot
        #---
        plt.close('all')
        #---
        #  Plot to visualize convergence of relative standard deviation (SD),
        #  skewness (SK) and kurtosis (KU) as a function of number of random
        #  samples
        #---
        SD = np.zeros_like(Powers2D[iso])
        SK = np.zeros_like(Powers2D[iso])
        KU = np.zeros_like(Powers2D[iso])
        nrand = np.shape(Powers2D[iso])[0]
        xaxis = np.arange(0, nrand)
        # Compute cumulative SD, SK and KU
        for i in xaxis:
            for x in range(0, len(CoreLayout[0, :])):
                for y in range(0, len(CoreLayout[:, 0])):
                    if CoreLayout[x, y] == 1:
                        SD[i, x, y] = (np.std(Powers2D[iso][:i + 1, x, y])
                                       /np.mean(Powers2D[iso][:i + 1, x, y])
                                       *100)
                        SK[i, x, y] = skew(Powers2D[iso][:i + 1, x, y])
                        KU[i, x, y] = kurtosis(Powers2D[iso][:i + 1, x, y])
        for stat in ['SD', 'SK', 'KU']:
            # Initialize figure
            fig, axs = plt.subplots(8, 8, sharex = 'all', sharey = 'all',
                                    figsize = (8, 8),
                                    gridspec_kw = {'hspace': 0, 'wspace': 0})
            for x in range(0, len(CoreLayout[0, :])):
                for y in range(0, len(CoreLayout[:, 0])):
                    if CoreLayout[x, y] == 0:
                        fig.delaxes(axs[x][y])
                    else:
                        if stat == 'SK' or stat == 'KU':
                            # Add a black vertical line on zero
                            axs[x, y].axhline(y = 0, linestyle = '--',
                                              linewidth = 0.5,
                                              color = 'k')
                        if stat == 'SD':
                            axs[x, y].plot(xaxis, SD[:, x, y])
                            # Zoom in the y axis on [0 ; 1.7%]
                            axs[x, y].set_ylim(0, 1.7)
                            # Add percent sign after yticks
                            labels = axs[x, y].get_yticks()
                            percentlabels = [0] + ['{:01.0f}'.format(i) + '%'
                                                   for i in labels[1:]]
                            axs[x, y].set_yticklabels(percentlabels)
                        elif stat == 'SK':
                            axs[x, y].plot(xaxis, SK[:, x, y])
                            # Zoom in the y axis
                            if iso.startswith('Ag') or iso == 'In115':
                                axs[x, y].set_ylim(-1.9, +1.9)
                            else:
                                axs[x, y].set_ylim(-0.65, +0.65)
                        elif stat == 'KU':
                            axs[x, y].plot(xaxis, KU[:, x, y])
                            # Zoom in the y axis
                            if iso == 'Ag107' or iso == 'In115':
                                axs[x, y].set_ylim(-3, +3)
                            else:
                                axs[x, y].set_ylim(-1.25, +1.25)
                        # Rotate labels
                        axs[x, y].tick_params(axis = 'x', labelrotation = 90)
            for x in range(0, len(CoreLayout[0, :])):
                for y in range(0, len(CoreLayout[:, 0])):
                    # Remove some of the axes, in order to avoid doubled axes
                    if y != len(CoreLayout[:, 0]) - 1:
                        # If there is a graph to its right, do not draw the
                        # right axis
                        if CoreLayout[x, y + 1] == 1:
                            axs[x, y].spines['right'].set_visible(False)
                    if x != len(CoreLayout[:, 0]) - 1:
                        # If there is a graph at his bottom, do not draw the
                        # bottom axis
                        if CoreLayout[x + 1, y] == 1:
                            axs[x, y].spines['bottom'].set_visible(False)
                    if x < 0:
                        axs[x, y].set_yticks([])
            #---
            #  Add a title
            #---
            if stat == 'SD':
                intitle = 'Relative standard deviation'
            elif stat == 'SK':
                intitle = 'Skewness'
            elif stat == 'KU':
                intitle = 'Excess kurtosis'
            st = fig.suptitle(intitle + ' of the assembly power as a function '
                              + 'of\nthe number of ' + iso + ' random samples '
                              + '(maximum ' + str(len(Powers2D[iso][:, 0, 0]))
                              + '), with Drakkar on\nTihange first '
                              + 'zero-power start-up, ' + text[controlrod])
            #---
            #  Save plot as pdf (vectorized)
            #---
            fig.savefig('output_2Dpow/' + stat + '_' + controlrod + '_' + iso
                        + '.pdf', bbox_extra_artists=[st], bbox_inches='tight')
            #---
            #  Clean-up for next plot
            #---
            plt.close('all')
    #---
    #  Merge pdf with ghostscript
    #---
    os.system('gs -q -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=output_2Dpow/'
              + controlrod + '.pdf -dBATCH output_2Dpow/' + controlrod
              + '_*.pdf')
#---
#  Plot, as a summary, the maximum relative standard deviation on an assembly
#  power, for each isotope and each control rod state
#---
w = 0.75
dimw = w/len(controlrods)
# Order according to std.dev. in the 'CD control rods inserted' state
Ordering = 'CD'
maxstd_sorted = {}
maxstd_sorted[Ordering] = sorted(maxstd[Ordering].items(),
                                 key = lambda x: x[1],
                                 reverse=True)
isotopes = [i[0] for i in maxstd_sorted[Ordering]]
relstd = [i[1] for i in maxstd_sorted[Ordering]]
# Sort in the same order (isotope-wise) the remaining states (except the one
# already sorted)
controlrods.remove(Ordering)
for CR in controlrods:
    maxstd_sorted[CR] = [maxstd[CR][i] for i in isotopes]
# Preparing a bar plot
fig, ax = plt.subplots()
y_pos = np.arange(len(isotopes))
# Plot the data
ax.barh(y_pos - 0.25, relstd, dimw, align='center',
        label = text[Ordering], zorder = 3)
i = 0
for CR in controlrods:
    i = i + 1
    ax.barh(y_pos + i * dimw - 0.25, maxstd_sorted[CR], dimw,
            align='center', label = text[CR].capitalize(), zorder = 3)
#---
#  Graph glitter
#---
ax.set_yticks(y_pos)
ax.set_yticklabels(isotopes)
ax.invert_yaxis() # labels read top-to-bottom
ax.set_ylim(len(isotopes) - 1 + w/2, -w/2)
ax.set_xlabel(r'Relative standard deviation (1$\sigma$) of the most '
              + 'uncertain assembly power')
plt.legend(loc = 'lower right', framealpha = 1.0)
ax.set_title('Ranking of nuclear data uncertainties with Drakkar on Tihange\n'
             + 'first zero-power start-up, with three control rods insertions')
# Add one minor tick between each major ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(n = 2))
# Add percent sign after yticks
labels = ax.get_xticks()
percentlabels = [0] + ['{:02.1f}'.format(i) + '%' for i in labels[1:]]
ax.set_xticklabels(percentlabels)
# Add a light x-axis grid in the background
ax.grid(which = 'major', axis = 'x', linewidth = 1.2, zorder = 0,
        linestyle = '--')
ax.grid(which = 'minor', axis = 'x', linewidth = 0.5, zorder = 0,
        linestyle = '--', dashes=(20, 20))
plt.tight_layout()
#---
#  Save plot as pdf (vectorized)
#---
fig.savefig('output_2Dpow/Summary.pdf')
#---
#  Clean-up
#---
plt.close('all')

print("Plotting completed")
