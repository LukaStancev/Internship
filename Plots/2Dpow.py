#
#  Plotting power uncertainty on 2D maps
#  Usage : python3 2Dpow.py
#

# Imports
import lcm
import numpy as np
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

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
    # Add the horizontal line. The information contained is redundant, except
    # for the central value
    LowerLine = np.hstack(( full[center, center], np.flip(Quarter[:, 0]) ))
    return np.vstack(( Quarter, LowerLine ))
# def unfold

#for controlrods in [ 'ARO', 'C', 'CD' ]:
os.system("ln -s ../Drakkar/Linux_x86_64/_PowerARO_*.ascii .")


# List isotopes we've been (potentially) randomly sampling
firstfile = glob.glob('_PowerARO_*.ascii')[0]
ResultFile = lcm.new('LCM_INP', firstfile[1:])
isotopes = ResultFile['NamIso'].split()

powers = {}
for iso in isotopes:
    powers[iso] = []

for file in list(glob.glob('_PowerARO_*.ascii')):
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
# Convert this list of 1D numpy array into a 2D numpy array: [irand, iassembly]
tmp = {}
for iso in isotopes:
    tmp[iso] = np.array(powers[iso])
powers = tmp
# Update 'isotopes' list so that only those isotopes that were indeed randomly
# sampled can be the subject of loops
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
#    print(Powers2D[iso])



#---
#  Plotting
#---


#print(powers.min())
#print(powers.max())
#print(len(powers))

### Exemple simple
#fig, ax = plt.subplots()
#ax.hist(powers[iso][:,0], bins=10)
#plt.savefig('myfig.png')
## Exemple 8x8
#fig, axs = plt.subplots(8,8)
#axs[0, 0].hist(powers[iso][:,0], bins=10)
#plt.savefig('myfig')
## Exemple 8x8 collé
#fig, axs = plt.subplots(8, 8, sharex='col', sharey='row',
#                        gridspec_kw={'hspace': 0, 'wspace': 0})
#axs[0, 0].hist(powers[iso][:,0], bins=10)
#plt.savefig('myfig')
# Exemple 8x8 collé sans les cases en trop

isotopes = ['U235']
#isotopes = ['U235', 'U238', 'O16']
for iso in isotopes:
    print(iso)
    fig, axs = plt.subplots(8, 8, sharex='all',
#                            constrained_layout=True,
#           subplot_kw=dict(adjustable='box'),
                            figsize=(8, 8),
                            gridspec_kw={'hspace': 0, 'wspace': 0})
#                                'height_ratios': [1]*8, 'width_ratios': [1]*8})
    for x in range(0, len(CoreLayout[0, :])):
        for y in range(0, len(CoreLayout[:, 0])):
            if CoreLayout[x, y] == 0:
                fig.delaxes(axs[x][y])
            else:
                # Plot an histogram with the assembly power
                axs[x, y].hist(Powers2D[iso][:, x, y], bins=20, density=True)
#                if x == 7 and (y == 0 or y == 7):
#                    mean = np.mean(Powers2D[iso][:, x, y])
#                    relstd = np.std(Powers2D[iso][:, x, y])/mean*100
#                #print(x, y, mean, relstd)
#                    print(iso, x, y, mean, relstd)


#                xlim = axs[x, y].get_xlim()
#                ylim = axs[x, y].get_ylim()
#                ratio = axs[x, y].get_data_ratio()
#                print(xlim, ylim, ratio)

                # Remove Y-ticks (a probability density does not require that)
                axs[x, y].set_yticks([])
    for x in range(0, len(CoreLayout[0, :])):
        for y in range(0, len(CoreLayout[:, 0])):
            # Assemblies are squares so we'd like square subplots
            axs[x, y].set_aspect(1./axs[x, y].get_data_ratio())

#    gs1 = GridSpec(8, 8, left=0.0, right=0.0)
#    for x in range(0, len(CoreLayout[0, :])):
#        for y in range(0, len(CoreLayout[:, 0])):
#            fig.update(gs1)

#    gs1 = gridspec.GridSpec(8, 8)
#    fig.update(wspace=0.0, hspace=0.) # set the spacing between axes.
#    fig.gridspec.GridSpec(hspace = 0.0, wspace = 0.0)

#    print('test...')
#    fig.subplots_adjust(wspace = 0.0, hspace = 0.0)
#    fig.tight_layout()
    fig.savefig(iso + '.pdf')

#    fig.savefig(iso + '.pdf', bbox_inches = 'tight', pad_inches = 0)


# Calculer moyenne/ecart-type et afficher ecart-type sur chaque case, sous la
# forme $\sigma=X.X%$

# Fond de couleur : moyenne de la case (min : min du min, max : max du max)
# Chaque baton d'une couleur diff = correspond a la puissance en question
# Largeur des distributions : max au centre et au bord, plus petit ailleurs
# Mais alors comment comparer les deux distributions...?

# Exemple en haut à droite, expliquant :
# Probability density vs power normalized on full core power

print("Plotting completed")
