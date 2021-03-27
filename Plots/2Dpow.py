#
#  Plotting power uncertainty on 2D maps
#  Usage : python3 2Dpow.py
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import lcm
import numpy as np
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from BasicFunctions import *
#
controlrods = ['CD', 'D', 'ARO']
text = {}
text['ARO'] = 'all rods out'
text['D'] = 'D rod bank inserted'
text['CD'] = 'C and D rod banks inserted'
# Create a dictionnary to store the maximum standard deviation among
# every assemblies
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
        fig, axs = plt.subplots(8, 8, sharex='all', figsize=(8, 8),
                                gridspec_kw={'hspace': 0, 'wspace': 0})
        for x in range(0, len(CoreLayout[0, :])):
            for y in range(0, len(CoreLayout[:, 0])):
                if CoreLayout[x, y] == 0:
                    fig.delaxes(axs[x][y])
                else:
                    # Plot an histogram with the assembly power
                    axs[x, y].hist(Powers2D[iso][:, x, y], bins=20,
                                   density=True)
                    # Compute relative standard deviation for each assembly
                    mean = np.mean(Powers2D[iso][:, x, y])
                    relstd = np.std(Powers2D[iso][:, x, y])/mean*100
                    # Store its maximum for future purposes
                    maxstd[controlrod][iso] = max(relstd,
                                                  maxstd[controlrod][iso])
                    # Print relative standard deviation for each assembly
                    textstr = r'$\sigma / \mu = %.1f$%%' % (relstd, )
                    axs[x, y].text(0.05, 0.95, textstr,
                                   transform=axs[x, y].transAxes,
                            #fontsize=10,
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
        #  Add a title
        #---
        fig.suptitle('Tihange first zero-power start-up, ' + text[controlrod]
                     + '\nProbability density vs normalized assembly power'
                     + '\n' + str(len(Powers2D[iso][:, 0, 0]))
                     + ' random samplings of ' + iso)
        #---
        #  Save plot as pdf (vectorized)
        #---
        os.system('mkdir -p output_2Dpow')
        fig.savefig('output_2Dpow/' + controlrod + '_' + iso + '.pdf')
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
ax.set_title('Ranking of nuclear data uncertainties on Tihange first\n'
             + 'zero-power start-up, with three control rods insertions')
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
