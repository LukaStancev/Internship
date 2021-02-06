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

def fold(full):
    # Retrieve center coordinates (assuming a square)
    center = len(full[:, 0])//2
    # We retrieve each quarter, that are not squares (+1 in one dimension)
    UpperRight = full[:center, center:]
    LowerRight = np.rot90(full[center:, center+1:])
    LowerLeft = np.rot90( full[center+1:, :center+1], 2)
    UpperLeft = np.rot90( full[:center+1, :center], 3)
    # Cumulate every quarters and mean them
    Quarter = (UpperRight + LowerRight + LowerLeft + UpperLeft)/4
    # Add the horizontal line. The information contained here is redundant,
    # except for the central value (which is full[center, center])
    LowerLine = np.hstack(( full[center, center], np.flip(Quarter[:, 0]) ))
    return np.vstack(( Quarter, LowerLine ))

for controlrods, controlrodstext in zip(['ARO', 'D', 'CD'], \
['all rods out', 'D rod bank inserted', 'C and D rod banks inserted']):
    print('controlrods = ' + controlrods)
    os.system('ln -s ../Drakkar/Linux_x86_64/_Power' + controlrods
              + '_*.ascii .')
    # List isotopes we've been (potentially) randomly sampling
    firstfile = glob.glob('_Power' + controlrods + '_*.ascii')[0]
    ResultFile = lcm.new('LCM_INP', firstfile[1:])
    isotopes = ResultFile['NamIso'].split()
    # Initialize dictionnary (isotope-wise) of lists (sampling-wise) that will
    # contain the retrieved data
    powers = {}
    for iso in isotopes:
        powers[iso] = []
    for file in list(glob.glob('_Power' + controlrods + '_*.ascii')):
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
    # Convert this list of 1D numpy array into a 2D numpy array:
    # [irand, iassembly]
    tmp = {}
    for iso in isotopes:
        tmp[iso] = np.array(powers[iso])
    powers = tmp
    # Update 'isotopes' list so that only those isotopes that were indeed
    # randomly sampled can be the subject of loops
    isotopes_used = []
    for iso in isotopes:
        if powers[iso].size != 0:
            isotopes_used.append(iso)
    isotopes = isotopes_used
    # Layout of the assemblies in the core:
    #   = 0, void (i.e. no assembly at this position)
    #   = 1, presence of an assembly
    if len(powers[iso][0, :]) == 157: # assemblies
        CoreLayout = np.array([[1, 1, 0, 0, 0, 0, 0, 0],
                               [1, 1, 1, 1, 0, 0, 0, 0],
                               [1, 1, 1, 1, 1, 0, 0, 0],
                               [1, 1, 1, 1, 1, 1, 0, 0],
                               [1, 1, 1, 1, 1, 1, 1, 0],
                               [1, 1, 1, 1, 1, 1, 1, 0],
                               [1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1]], dtype=int)
    else:
        raise Exception('Unsupported core layout')
    # --- Unfold this quarter core on a full core layout
    # Prepare upper core layout, minus central vertical line
    UpperLeft = np.delete( np.rot90(CoreLayout), -1, axis=0)
    UpperRight = np.delete( np.delete(CoreLayout, -1, axis=0), 0, axis=1 )
    Upper = np.hstack( (UpperLeft, UpperRight) )
    # Prepare lower core layout, letting vertical line
    LowerLeft = np.rot90(CoreLayout, 2)
    LowerRight = np.delete( np.rot90(CoreLayout, 3), 0, axis=1 )
    Lower = np.hstack( (LowerLeft, LowerRight) )
    FullCoreLayout = np.vstack((Upper, Lower))
    #
    FullPowers2D = {}
    Powers2D = {}
    for iso in isotopes:
        # Deploy the 1D-indexed power maps into a 2D-indexed power maps
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
        fig.suptitle('Tihange first zero-power start-up, ' + controlrodstext
                     + '\nProbability density vs normalized assembly power'
                     + '\n' + str(len(Powers2D[iso][:, 0, 0]))
                     + ' random samplings of ' + iso)
        #---
        #  Save plot as pdf (vectorized)
        #---
        fig.savefig(controlrods + '_' + iso + '.eps')
        #---
        #  Clean-up for next plot
        #---
        plt.close('all')
    #---
    #  Merge pdf with ghostscript
    #---
    os.system('gs -q -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=' + controlrods
              + '.pdf -dBATCH ' + controlrods + '_*.pdf')
print("Plotting completed")
