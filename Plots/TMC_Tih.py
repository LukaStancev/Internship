#
#  Plotting power on 2D maps for TMC samples, on Tihange
#  Usage : python3 TMC_Tih.py
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import lcm
import numpy as np
from scipy.stats import skew
from scipy.stats import kurtosis
import math
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.gridspec import GridSpec
from BasicFunctions import *
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en

# Initialize a gaussian vector (used for legends)
gaussian = np.random.randn(300)
# Declare control rod insertions and its full names
controlrods = ['CD', 'D', 'ARO']
text = {}
if lang == 'en':
    text['ARO'] = 'all rods out'
    text['D'] = 'D rod bank inserted'
    text['CD'] = 'C and D rod banks inserted'
elif lang == 'fr':
    text['ARO'] = 'toutes grappes extraites'
    text['D'] = 'groupe D inséré'
    text['CD'] = 'groupes C et D insérés'
# Create a dictionary to store the standard deviation of the assembly with the
# most uncertain power, among all those assemblies in the core.
maxstd = {}
# ...and it's standard error
maxstderr = {}
# Same for maximum/minimum skewness and kurtosis
minSK = {}
maxSK = {}
minKU = {}
maxKU = {}
nsamples = {}
for controlrod in controlrods:
    print('controlrod = ' + controlrod)
    maxstd[controlrod] = {}
    maxstderr[controlrod] = {}
    minSK[controlrod] = {}
    maxSK[controlrod] = {}
    minKU[controlrod] = {}
    maxKU[controlrod] = {}
    # Create links toward every available power distribution
    os.system('ln -s ../Drakkar/Output_TIH_TMC/_Power' + controlrod
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
            if int(ResultFile[iso]) != -33:
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
        nsamples[iso] = len(Powers2D[iso][:, 0, 0])
        # Initialize maximum relative standard deviation to zero
        maxstd[controlrod][iso] = 0
        maxstderr[controlrod][iso] = 0
        minSK[controlrod][iso] = math.inf
        maxSK[controlrod][iso] = -math.inf
        minKU[controlrod][iso] = math.inf
        maxKU[controlrod][iso] = -math.inf
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
                    # Plot an histogram with the assembly power
                    axs[x, y].hist(Powers2D[iso][:, x, y], bins = 20,
                                   density = True)
                    # Compute relative standard deviation for each assembly
                    mean = np.mean(Powers2D[iso][:, x, y])
                    relstd = np.std(Powers2D[iso][:, x, y], ddof = 1)/mean*100
                    textstr = (r'$\mu$=$%.2f$' % mean + '\n'
                               + r'$\sigma$=$%.1f$' % (relstd, ) + '\%')
                    # Store its maximum for future purposes (summary plot)
                    if relstd > maxstd[controlrod][iso]:
                        maxstd[controlrod][iso] = relstd
                        SE = SE_SD(Powers2D[iso][:, x, y])
                        maxstderr[controlrod][iso] = SE/mean*100
                    # Compute skewness and kurtosis only if sigma is large
                    # enough (> 0.1%)
                    if relstd > 0.1:
                        skwn = skew(Powers2D[iso][:, x, y])
                        kurt = kurtosis(Powers2D[iso][:, x, y])
                        textstr += '\n$S$=$%.1f$' % skwn
                        textstr += '\n$K$=$%.1f$' % kurt
                        if skwn < minSK[controlrod][iso]:
                            minSK[controlrod][iso] = skwn
                        if skwn > maxSK[controlrod][iso]:
                            maxSK[controlrod][iso] = skwn
                        if kurt < minKU[controlrod][iso]:
                            minKU[controlrod][iso] = kurt
                        if kurt > maxKU[controlrod][iso]:
                            maxKU[controlrod][iso] = kurt
                    # Print statistics for each assembly
                    axs[x, y].text(0.05, 0.95, textstr,
                                   transform=axs[x, y].transAxes,
                                   verticalalignment = 'top')
                    # Remove Y-ticks (unnecessary for a probability density)
                    axs[x, y].set_yticks([])
        # If no skewness or kurtosis has been calculated (due to insufficient
        # standard deviation), then it makes no sense to plot them. In this
        # case, we set them to zero.
        if minSK[controlrod][iso] == math.inf:
            minSK[controlrod][iso] = 0
            maxSK[controlrod][iso] = 0
            minKU[controlrod][iso] = 0
            maxKU[controlrod][iso] = 0
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
            textstr = (r'$\mu$ : mean of the power assembly' + '\n'
                       + r'$\sigma$ : relative standard deviation '
                       + '(in \%)\n'
                       + r'$S$ : skewness (shown if $\sigma > 0.1$'
                       + '\%)\n'
                       + r'$K$ : excess kurtosis (if $\sigma > 0.1$'
                       + '\%)')
        elif lang == 'fr':
            textstr = (r'$\mu$ : puissance moyenne de ' + 'l\'assemblage\n'
                       + r'$\sigma$ : écart-type relatif '
                       + '(en \%)\n'
                       + r'$S$ : asymétrie (affichée si $\sigma > 0.1$'
                       + '\%)\n'
                       + r'$K$ : excès de kurtosis (si $\sigma > 0.1$'
                       + '\%)')
        axs[0, 2].text(0.05, 1.04, textstr, verticalalignment = 'top',
                       transform = axs[0, 2].transAxes)
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
        #  Add a title and save plot as pdf (vectorized)
        #---
        if lang == 'en':
            title = ('Probability density vs. normalized assembly power, with '
                     + 'Drakkar and\n'
                     + 'JEFF-3.3 best estimate except for '
                     + iso.replace('_', '\_') + ' ('
                     + str(len(Powers2D[iso][:, 0, 0]))
                     + ' random samples) on\nTihange first zero-power'
                     + ' start-up, ' + text[controlrod])
        elif lang == 'fr':
            title = ''
        os.system('mkdir -p output_TMC_Tih')
        if title:
            st = fig.suptitle(title)
            fig.savefig('output_TMC_Tih/' + controlrod + '_' + iso + '.pdf',
                        bbox_extra_artists=[st], bbox_inches='tight')
        else:
            fig.savefig('output_TMC_Tih/' + controlrod + '_' + iso + '.pdf',
                        bbox_inches='tight')
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
                        SD[i, x, y] = (np.std(Powers2D[iso][:i + 1, x, y],
                                       ddof = 1)
                                       /np.mean(Powers2D[iso][:i + 1, x, y])
                                       *100)
                        SK[i, x, y] = skew(Powers2D[iso][:i + 1, x, y])
                        KU[i, x, y] = kurtosis(Powers2D[iso][:i + 1, x, y])
        for stat in ['SD', 'SK', 'KU']:
            # Initialize figure
            fig, axs = plt.subplots(8, 8, sharex = 'all', sharey = 'all',
                                    figsize = set_size('square', bonus = True),
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
                            percentlabels = [0] + ['{:01.0f}'.format(i) + '\%'
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
            #  Add a title and save plot as pdf (vectorized)
            #---
            if lang == 'en':
                if stat == 'SD':
                    intitle = 'Relative standard deviation'
                elif stat == 'SK':
                    intitle = 'Skewness'
                elif stat == 'KU':
                    intitle = 'Excess kurtosis'
                title = (intitle + ' of the assembly power as a function '
                         + 'of\nthe number of ' + iso.replace('_', '\_')
                         + ' random samples (maximum '
                         + str(len(Powers2D[iso][:, 0, 0]))
                         + '), with Drakkar on\nTihange first '
                         + 'zero-power start-up, ' + text[controlrod])
            elif lang == 'fr':
                title = ''
            if title:
                st = fig.suptitle(title)
                fig.savefig('output_TMC_Tih/' + stat + '_' + controlrod + '_'
                            + iso + '.pdf', bbox_inches='tight',
                            bbox_extra_artists=[st])
            else:
                fig.savefig('output_TMC_Tih/' + stat + '_' + controlrod + '_'
                            + iso + '.pdf', bbox_inches='tight')
            #---
            #  Clean-up for next plot
            #---
            plt.close('all')
    #---
    #  Merge pdf with ghostscript
    #---
    os.system('gs -q -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=output_TMC_Tih/'
              + controlrod + '.pdf -dBATCH output_TMC_Tih/' + controlrod
              + '_*.pdf')
#---
#  For each isotope and each control rod state, plot, as a summary:
#  * the maximum relative standard deviation of an assembly power
#  * the minimum and maximum skewness on the core
#  * the minimum and maximum kurtosis on the core
#  The standard errors (due to the finite statistics) are plotted for each of
#  these.
#---
#---
#  Sorting standard deviations
#---
# Order according to std.dev. in the 'CD control rods inserted' state
Ordering = 'CD'
maxstd_sorted = {}
maxstd_sorted[Ordering] = sorted(maxstd[Ordering].items(),
                                 key = lambda x: x[1],
                                 reverse=True)
isotopes = [i[0] for i in maxstd_sorted[Ordering]]
relstd = [i[1] for i in maxstd_sorted[Ordering]]
# Sort in the same order (isotope-wise) the std.dev. of the other control rod
# states, and also the standard errors for all control rod states
maxstderr_sorted = {}
controlrods.remove(Ordering)
for CR in controlrods:
    maxstd_sorted[CR] = [maxstd[CR][i] for i in isotopes]
    maxstderr_sorted[CR] = [maxstderr[CR][i] for i in isotopes]
maxstderr_sorted[Ordering] = [maxstderr[Ordering][i] for i in isotopes]
#---
#  Graphical parameters
#---
w = 0.75
dimw = w/(len(controlrods) + 1)
# Proportionality factor for the cap sizes of the error bars
nbar = (len(controlrods) + 2)*len(isotopes) - 1
capsize = 125/nbar
#---
#  Create subplots with proper sizes
#---
# Filter the isotopes so that the skewness and kurtosis are not displayed when
# the relative standard deviation is smaller than 0.1%.
isotopes_maj = isotopes[0:np.count_nonzero(np.array(relstd) > 0.1)]
nbar_maj = (len(controlrods) + 2)*len(isotopes_maj) - 1
fig = plt.figure(figsize = set_size())
gs = GridSpec(2, 2, figure = fig, wspace = 0.0, hspace = 0.0,
              height_ratios = [nbar_maj, nbar - nbar_maj])
#---
#  Plot standard deviations
#---
ax1 = fig.add_subplot(gs[:, 0])
y_pos = np.arange(len(isotopes))
label = text[Ordering][0].capitalize() + text[Ordering][1:]
ax1.barh(y_pos - 0.25, relstd, dimw,
         xerr = maxstderr_sorted[Ordering],
         align = 'center', label = label, zorder = 3,
         error_kw = dict(capsize = capsize, capthick = w*2/3, lw = w*2/3))
i = 0
for CR in controlrods:
    i = i + 1
    label = text[CR][0].capitalize() + text[CR][1:]
    ax1.barh(y_pos + i * dimw - 0.25, maxstd_sorted[CR], dimw,
             xerr = maxstderr_sorted[CR],
             align = 'center', label = label, zorder = 3,
             error_kw = dict(capsize = capsize, capthick = w*2/3, lw = w*2/3))
# Reset controlrods list
controlrods = [Ordering] + controlrods
#---
#  Graph glitter
#---
ax1.set_yticks(y_pos)
ax1.set_yticklabels([iso.replace('_', '\_') for iso in isotopes])
ax1.invert_yaxis() # labels read top-to-bottom
ax1.set_ylim(len(isotopes) - 1 + w/2, -w/2)
if lang == 'en':
    ax1.set_xlabel(r'Relative standard deviation (1$\sigma$) of'
                   + '\nthe most uncertain assembly power')
elif lang == 'fr':
    ax1.set_xlabel(r'Écart-type relatif (1$\sigma$) de la puissance'
                   + '\nd\'assemblage la plus incertaine')
plt.legend(loc = 'lower center', framealpha = 1.0)
# Add one minor tick between each major ticks
ax1.xaxis.set_minor_locator(AutoMinorLocator(n = 2))
# Add percent sign after yticks
labels = ax1.get_xticks()
percentlabels = [0] + ['{:02.1f}'.format(i) + '\%' for i in labels[1:]]
ax1.set_xticklabels(percentlabels)
# Add a light x-axis grid in the background
ax1.grid(which = 'major', axis = 'x', linewidth = 1.2, zorder = 0,
        linestyle = '--')
ax1.grid(which = 'minor', axis = 'x', linewidth = 0.5, zorder = 0,
        linestyle = '--', dashes=(20, 20))
#---
#  Skewness and kurtosis
#---
# Add a common title below skewness and kurtosis (produced with a nested
# gridspec), hiding a phantom figure
axtitle = fig.add_subplot(gs[0, -1])
axtitle.set_yticklabels([])
axtitle.set_yticks([])
axtitle.spines['top'].set_visible(False)
axtitle.spines['bottom'].set_visible(False)
axtitle.spines['left'].set_visible(False)
axtitle.spines['right'].set_visible(False)
axtitle.set_xticks([0])
axtitle.set_xticklabels(['0'])
axtitle.tick_params(axis = 'x', colors = 'white')
axtitle.set_xlim([-1, 1])
if lang == 'en':
    axtitle.set_xlabel('Higher moments of\n'
                       + 'probability distributions:\n'
                       + 'minimum/maximum\n'
                       + 'skewness (left) and excess\n'
                       + 'kurtosis (right) of assembly power\n'
                       + r'(displayed only if $\sigma > 0.1$' + '\%)')
elif lang == 'fr':
    axtitle.set_xlabel('Moments supérieurs des\n'
                       + 'distributions de probabilité :\n'
                       + 'minimum/maximum de\n'
                       + 'l\'asymétrie (à gauche) et de\n'
                       + 'l\'excès de kurtosis (à droite)\n'
                       + 'de la puissance par assemblage\n'
                       + r'(affichés si $\sigma > 0.1$' + '\%)')
# For each isotope, compute standard errors of skewness and kurtosis
SE_skwn = [SE_SK(nsamples[iso]) for iso in isotopes_maj]
SE_kurt = [SE_KU(nsamples[iso]) for iso in isotopes_maj]
# Duplicate the standard errors, in order to apply them to both the min and max
SE_skwn = SE_skwn + SE_skwn
SE_kurt = SE_kurt + SE_kurt
# Prepare the subplots with a nested gridspec
inner_grid = gs[0, -1].subgridspec(1, 2, wspace = 0, hspace = 0)
ax2 = fig.add_subplot(inner_grid[0, -2])
ax3 = fig.add_subplot(inner_grid[0, -1])
for jCR, CR in enumerate(controlrods):
    # Align the points horizontally with the standard deviation barplots
    height = -(np.arange(0, len(isotopes_maj))*(len(controlrods) + 1)
               + 1/2 + jCR)/nbar_maj
    # Prepare the data to be plotted (min and max)
    skwn = ([minSK[CR][iso] for iso in isotopes_maj]
            + [maxSK[CR][iso] for iso in isotopes_maj])
    kurt = ([minKU[CR][iso] for iso in isotopes_maj]
            + [maxKU[CR][iso] for iso in isotopes_maj])
    heights = list(height) + list(height)
    # Filter out values equal to zero (unattributed) that can happen in 'all
    # rods out' state with isotopes present only in control rods
    nonzero = np.invert(np.isclose(skwn, 0))
    heights = np.array(heights)[nonzero]
    skwn = np.array(skwn)[nonzero]
    kurt = np.array(kurt)[nonzero]
    # Plot all of this data
    ax2.errorbar(skwn, heights, fmt = 'x', xerr = np.array(SE_skwn)[nonzero],
                 capsize = capsize, capthick = w*2/3, elinewidth = w*2/3,
                 ms = 2)
    ax3.errorbar(kurt, heights, fmt = 'x', xerr = np.array(SE_kurt)[nonzero],
                 capsize = capsize, capthick = w*2/3, elinewidth = w*2/3,
                 ms = 2)
#---
#  Graph glitter
#---
# Remove useless labels and ticks on height
ax2.set_yticklabels([])
ax2.set_yticks([])
# Add isotopes labels, on the right
height = -(np.arange(0, len(isotopes_maj))*(len(controlrods) + 1)
           + 1/2 + int((len(controlrods)-1)/2))/nbar_maj
ax3.yaxis.set_ticks_position('right')
ax3.set_yticks(height)
ax3.set_yticklabels([iso.replace('_', '\_') for iso in isotopes_maj])
# Show proper y limits to see the points in the right positions
ax2.set_ylim([-1, 0])
ax3.set_ylim([-1, 0])
# Set kurtosis major ticks every 2
ax3.xaxis.set_major_locator(MultipleLocator(2))
# Add one minor tick between each major ticks
ax2.xaxis.set_minor_locator(AutoMinorLocator(n = 2))
ax3.xaxis.set_minor_locator(AutoMinorLocator(n = 2))
# And labels on minor ticks of kurtosis, with ticks of the same size
ax3.tick_params(axis = 'x', which = 'minor',
                length = plt.rcParams["xtick.major.size"])
ax3.xaxis.set_minor_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
# Add a light x-axis grid in the background
ax2.grid(which = 'major', axis = 'x', linewidth = 1.2, zorder = 0,
        linestyle = '--')
ax2.grid(which = 'minor', axis = 'x', linewidth = 0.5, zorder = 0,
        linestyle = '--', dashes=(20, 20))
ax3.grid(which = 'major', axis = 'x', linewidth = 1.2, zorder = 0,
        linestyle = '--')
ax3.grid(which = 'minor', axis = 'x', linewidth = 0.5, zorder = 0,
        linestyle = '--', dashes=(20, 20))
#---
#  Add a global title and save plot as pdf (vectorized)
#---
if lang == 'en':
    title = ('Ranking of nuclear data uncertainties with Drakkar on Tihange\n'
             + 'first zero-power start-up, with three control rods insertion\n'
             + 'All error bars correspond to statistical standard errors '
             + r'(1$\sigma$)')
elif lang == 'fr':
    title = ''
if title :
    st = fig.suptitle(title, y = 1.01)
    fig.savefig('output_TMC_Tih/Summary.pdf', bbox_inches = 'tight',
                extra_artists = [st])
else:
    fig.savefig('output_TMC_Tih/Summary.pdf', bbox_inches = 'tight')
#---
#  Clean-up
#---
plt.close('all')
print('Rough magnitude order of the maximum total uncertainty:')
for CR in controlrods:
    ssum = np.sum(list(maxstd[CR].values()))
    print(text[CR] + ', 3 sigma, arithmetic sum = ' + str(3*ssum))
    qsum = math.sqrt(np.sum(np.array(list(maxstd[CR].values()))**2))
    print(text[CR] + ', 2 sigma, quadratic sum = ' + str(2*qsum))
print("Plotting completed")
