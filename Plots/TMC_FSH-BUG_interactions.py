#
#  Plotting detector relative difference on 2D maps
#  Usage : python3 TMC_FSH-BUG_interactions.py
#  Author ; V. Salino (IRSN), 09/2021
#

# Imports
import lcm
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import skew
from scipy.stats import kurtosis
import os
import math
import subprocess
import glob
import re
import csv
import string
import pickle
import itertools
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
from BasicFunctions import *
import serpentTools
from serpentTools.settings import rc
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en

#---
#  Initialize the container of all combinations of detector responses
#---
responses = {}
# Layout of the assemblies in the core
CoreLayout, FullCoreLayout = GetCoreLayout(157)
FullCoreLayout4X = np.repeat(np.repeat(FullCoreLayout,
                                       2, axis = 0), 2, axis = 1)
#---
#  Retrieve responses measured experimentally by 50 mobile fission chambers
#---
for file in (set(glob.glob('../Measurements/CP0/*'))
# Minus Almaraz, which contains only measurements-derived power distribution
             - set(glob.glob('../Measurements/CP0/Almaraz*'))):
    print('Reading ' + file)
    lines = []
    with open(file) as csvfile:
        reader = csv.reader(decomment(csvfile))
        for line in reader:
            lines.append(line)
    # Replace spaces with NaNs
    respon = pd.DataFrame(lines).replace(r'^\s*$', np.nan, regex=True)
    # Transform values from strings to floats
    respon = respon.astype('float')
    # Check that normalization is properly done, up to a criterion (0.05%)
    nelem = np.sum(respon.count())
    totsum = respon.sum().sum()
    if np.abs(totsum/nelem - 1) > 0.05/100:
        raise Exception('The normalization of the file ' + file
                        + ' is incorrect. The sum of the responses is '
                        + str(totsum) + ' instead of ' + str(nelem)
                        + '. Review its responses.')
    # Renormalize to remove any possible small remaining residue
    respon = respon*nelem/totsum
    # Identify the source
    reactor = re.split('[/.]', file)[-2]
    responses[reactor] = np.array(respon)
#---
#  Include the average of the experimental measurements
#---
totsum = 0
for reactor in responses:
    totsum = totsum + responses[reactor]
responses['ExpMean'] = totsum/len(responses)
#---
#  Retrieve best-estimate Serpent power distribution
#---
# Shorter execution:
links = glob.glob('../Serpent/SuperBatch/Fessenheim-Bugey*4.sss2_det0.m')
# Retrieve Serpent power distribution
links = glob.glob('../Serpent/SuperBatch/Fessenheim-Bugey*.sss2_det0.m')
# Storage for independant runs' results
indruns = []
for link in links:
    indruns.append(ReadPdistrSerpent(link, FullCoreLayout))
powerSerpent = Deploy2D(np.mean(indruns, axis = 0), FullCoreLayout)
#---
#  Retrieve best-estimate Drakkar power distributions
#---
path = '../Drakkar/Output_FSH-BUG_BestEstimate/'
link = path + '_PowerARO_7bd89e1fe9cdb4e037a9ca84536a142d.ascii'
power = loadlcm(link)['POWER-CHAN']
# Each FullCoreLayout value is duplicated in 4 values
power = matrixaverage(Deploy2D(power, FullCoreLayout4X), 2)
# Filter out reflectors (=0 in FullCoreLayout) with NaNs
power = np.where(FullCoreLayout == 0, np.nan, power)
# Normalize
powerDrakkar = power/np.nanmean(power)
#---
#  Compute prior khi square
#---
respon = Deploy2D(loadlcm(link)['RESPON'], FullCoreLayout)
respon = respon*powerSerpent/powerDrakkar
# Filter out unmeasured responses with NaNs, in order to have a
# normalization identical to the measurements, i.e.
# sum(50 measured positions)=50.
respon_filtered = np.where(np.isnan(responses['Fessenheim-2']), np.nan,
                          respon)
normfactor = np.nanmean(respon_filtered)
# Filter out the reflectors (=0 in FullCoreLayout) with NaNs and apply that
# normalization
respon = np.where(FullCoreLayout == 0, np.nan, respon)/normfactor
# Quick overview of Khi2 contributions, with best estimate nuclear data
std = 0.02
perassembly = ((respon - responses['ExpMean'])/std)**2
khisqrs = np.nansum(perassembly)
print("Khi2 :")
print(khisqrs)
print("Khi2 contributions :")
print(pd.DataFrame(perassembly))
#---
#  Retrieve randomly sampled Drakkar powers and responses
#---
# Shorter execution:
links = glob.glob("../Drakkar/Output_FSH-BUG_TMC_interactions/_PowerARO_1*")
links = glob.glob("../Drakkar/Output_FSH-BUG_TMC_interactions/_PowerARO_*")
responsesTMC = []
responsesfTMC = []
powersTMC = []
powersuTMC = []
# Load any file, just to get the sampled isotopes
lcmfile = loadlcm(links[0])
isotopes = lcmfile['NamIso'].split()
indices = {}
for iso in isotopes:
    index = lcmfile[iso]
    if index != -33:
        indices[iso] = []
for link in links:
    lcmfile = loadlcm(link)
    # Retrieve power
    power = matrixaverage(Deploy2D(lcmfile['POWER-CHAN'], FullCoreLayout4X), 2)
    power = np.where(FullCoreLayout == 0, np.nan, power)
    # Normalize
    power = power/np.nanmean(power)
    # Correct for the deterministic bias, using Serpent results
    power = power*powerSerpent/powerDrakkar
    powersTMC.append(fold(power))
    powersuTMC.append(power)
    # Retrieve and deploy in 2D the Drakkar responses
    respon = Deploy2D(loadlcm(link)['RESPON'], FullCoreLayout)
    # Correct for the deterministic bias, using Serpent results
    respon = respon*powerSerpent/powerDrakkar
    # Filter out unmeasured responses with NaNs, in order to have a
    # normalization identical to the measurements, i.e.
    # sum(50 measured positions)=50.
    respon_filtered = np.where(np.isnan(responses['Fessenheim-2']), np.nan,
                              respon)
    normfactor = np.nanmean(respon_filtered)
    # Filter out the reflectors (=0 in FullCoreLayout) with NaNs and apply that
    # normalization
    respon = np.where(FullCoreLayout == 0, np.nan, respon)/normfactor
    responsesTMC.append(respon)
    responsesfTMC.append(fold(respon))
    # For each (sampled) isotope, retrieve the sample number
    isotopes = lcmfile['NamIso'].split()
    for iso in isotopes:
        index = int(lcmfile[iso])
        # Perform an action only if that isotope has been sampled
        if index != -33:
            indices[iso].append(index)
# Convert this list of 1D numpy array into a 2D numpy array
responsesTMC = np.array(responsesTMC)
responsesfTMC = np.array(responsesfTMC)
powersTMC = np.array(powersTMC)
powersuTMC = np.array(powersuTMC)
#---
# Based on
# * Rochman, D., Bauge, E., Vasiliev, A. et al. "Monte Carlo nuclear data
#   adjustment via integral information." Eur. Phys. J. Plus 133, 537 (2018).
#   https://tendl.web.psi.ch/bib_rochman/epjp1800496-offprints.pdf
# * P. Helgesson, H. Sjöstrand, A.J. Koning, J. Rydén, D. Rochman,
#   E. Alhassan and S. Pomp, "Combining Total Monte Carlo and Unified Monte
#   Carlo: Bayesian nuclear data uncertainty quantification from
#   auto-generated experimental covariances." Progress in Nuclear Energy,
#   Volume 96 (2017).
#   https://doi.org/10.1016/j.pnucene.2016.11.006
# * Rochman, D., Vasiliev, A., Ferroukhi, H. et al. "Correlation nu_p-sigma for
#   U-Pu in the thermal and resonance neutron range via integral information."
#   Eur. Phys. J. Plus 134, 453 (2019).
#   https://tendl.web.psi.ch/bib_rochman/epjp1900495-offprints.pdf
#---
nmes = np.count_nonzero(~np.isnan(responses['ExpMean']))
ncalc = len(responsesTMC)
expuncertainty = 0.02 # 2% (expert opinion)
errperassembly = np.copy(responses['ExpMean'])
errperassembly[~np.isnan(responses['ExpMean'])] = 0
khisqrs = np.empty(ncalc)
i = 0
# Compute khi-squares for each nuclear data sample
for i in range(0, ncalc):
    elm = ((responsesTMC[i] - responses['ExpMean'])/expuncertainty)**2
    khisqrs[i] = np.nansum(elm)
    errperassembly = errperassembly + elm
print("Khi2 :")
print(stats.describe(khisqrs))
print("Khi2 contributions :")
print(pd.DataFrame(errperassembly))
for nbins in [100, 150, 200, 250, 300, 350, 400]:
    fig, ax = plt.subplots()
    hist, bins = np.histogram(khisqrs, bins=nbins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    ax.hist(khisqrs, bins=logbins)
    plt.xscale('log')
    ax.set_xlabel(r'$\chi^2$')
    ax.set_ylabel('Nombre d\'échantillons par colonne')
    if lang == 'en':
        ax.set_title(str(nbins) + r'-bins histogram of detectors $\chi^2$')
    # Add various statistical informations
    stat = khisqrs
    text = '\n'.join(('Minimum=' + '{:.2f}'.format(np.min(stat)),
                      'Maximum=' + '{:.2f}'.format(np.max(stat)),
                      'Moyenne=' + '{:.2f}'.format(np.mean(stat)),
                      'Médiane=' + '{:.2f}'.format(np.quantile(stat, 0.5)),
                      'Écart-type=' + '{:.2f}'.format(np.std(stat)),
                      'Asymétrie=' + '{:.2f}'.format(skew(stat)),
                      'Excès de kurtosis=' + '{:.2f}'.format(kurtosis(stat))))
    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.65, 0.95, text, transform=ax.transAxes,
                    verticalalignment='top', bbox=props)
    fig.savefig('output_TMC_FSH-BUG_interactions/Khi2' + '_' + str(nbins)
                + 'bins.pdf', bbox_inches='tight')
    # Clean-up for next plot
    plt.close('all')
    del fig
    del ax
for method in ['Prior', 'BMC', 'BFMC']:
    if method == 'Prior':
        weights = np.ones_like(khisqrs)
    elif method == 'BMC':
        weights = np.exp(-khisqrs/2)
    elif method == 'BFMC':
        weights = np.exp(-khisqrs/np.min(khisqrs))
        weights_BFMC = weights
    # Normalize weights
    weights = weights/sum(weights)
    print(method + ' weights :')
    print(stats.describe(weights))
    # Save weights for later use
    if method == 'BMC':
        weights_BMC = weights
    elif method == 'BFMC':
        weights_BFMC = weights
    for nbins in [100, 150, 200, 250, 300, 350, 400]:
        fig, ax = plt.subplots()
        hist, bins = np.histogram(weights, bins=nbins)
        logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        ax.hist(weights, bins=logbins)
        plt.xscale('log')
        ax.set_xlabel(r'Poids $\omega_i$')
        ax.set_ylabel('Nombre d\'échantillons par colonne')
        if lang == 'en':
            ax.set_title(str(nbins) + r'-bins histogram of ' + method
                         + ' weights')
        # Add various statistical informations
        stat = weights
        text = '\n'.join(('Minimum=' + '{:.2e}'.format(np.min(stat)),
                          'Maximum=' + '{:.2e}'.format(np.max(stat)),
                          'Moyenne=' + '{:.2e}'.format(np.mean(stat)),
                          'Médiane=' + '{:.2e}'.format(np.quantile(stat, 0.5)),
                          'Écart-type=' + '{:.2e}'.format(np.std(stat)),
                          'Asymétrie=' + '{:.2f}'.format(skew(stat)),
                          'Excès de kurtosis=' + '{:.2f}'.format(kurtosis(stat)),
                          # ESS = Equivalent Sample Size
                          'ESS=' + '{:.2f}'.format(np.sum(weights)**2/
                                                   np.sum(weights**2))))
        # Some BMC weights can be very low, but that's not so interesting
        ax.set_xlim(10**(-10), 10**(-2))
        # place a text box in upper left in axes coords
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, text, transform=ax.transAxes,
                        verticalalignment='top', bbox=props)
        fig.savefig('output_TMC_FSH-BUG_interactions/Weights' + '_' + method
                    + '_' + str(nbins) + 'bins.pdf', bbox_inches='tight')
        # Clean-up for next plot
        plt.close('all')
        del fig
        del ax
    # Prior/posterior power average
    print(pd.DataFrame(np.average(powersTMC, axis=0)/
                       np.average(powersTMC, axis=0, weights = weights)-1)*100)
    #---
    #  Plotting prior and posterior power distributions and detector reponses
    #---
    for qty in ['pow', 'det']:
        # Initialize figure
        fig, axs = plt.subplots(len(CoreLayout), len(CoreLayout),
                                sharex = 'all',
                                figsize = set_size('square'),
                                gridspec_kw = {'hspace': 0, 'wspace': 0})
        for x in range(0, len(CoreLayout[0, :])):
            for y in range(0, len(CoreLayout[:, 0])):
                if CoreLayout[x, y] == 0:
                    axs[x, y].set_axis_off()
                else:
                    # Select quantity to be plotted
                    if qty == 'pow':
                        distrib = powersTMC[:, x, y]
                    elif qty == 'det':
                        distrib = responsesfTMC[:, x, y]
                    # Plot an histogram with the assembly power
                    axs[x, y].hist(distrib, bins = 20, density = True,
                                   weights = weights)
                    # Compute relative standard deviation for each assembly
                    mean = np.average(distrib, weights = weights)
                    var = np.average((distrib - mean)**2,
                                     weights = weights)
                    relstd = np.sqrt(var)/mean*100
                    textstr = (r'$\mu$=$%.2f$' % mean + '\n'
                               + r'$\sigma$=$%.1f$' % (relstd, ) + '\%')
                    # Compute skewness and kurtosis only if sigma is large
                    # enough (> 0.1%)
                    if relstd > 0.1:
                        skwn = weighted_skew(distrib, weights)
                        kurt = weighted_kurtosis(distrib, weights)
                        textstr += '\n$S$=$%.1f$' % skwn
                        textstr += '\n$K$=$%.1f$' % kurt
                    # Print statistics for each assembly
                    axs[x, y].text(0.05, 0.95, textstr, fontsize = 12,
                                   transform=axs[x, y].transAxes,
                                   verticalalignment = 'top')
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
            if qty == 'pow':
                textstr = r'$\mu$ : mean of the power assembly' + '\n'
            elif qty == 'det':
                textstr = r'$\mu$ : mean of the detector response' + '\n'
            textstr += (r'$\sigma$ : relative standard deviation '
                        + '(in \%)\n'
                        + r'$S$ : skewness (shown if $\sigma > 0.1$'
                        + '\%)\n'
                        + r'$K$ : excess kurtosis (if $\sigma > 0.1$'
                        + '\%)')
        elif lang == 'fr':
            if qty == 'pow':
                textstr = (r'$\mu$ : puissance moyenne de '
                           + 'l\'assemblage\n')
            elif qty == 'det':
                textstr = (r'$\mu$ : réponse moyenne du détecteur' + '\n')
            textstr += (r'$\sigma$ : écart-type relatif '
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
        gaussian = np.random.randn(300)
        axs[1, 6].set_axis_on()
        axs[1, 6].hist(1.0225 + 0.015 * gaussian, bins = 20, density = True)
        axs[1, 6].set_yticks([])
        axs[1, 6].set_aspect(1./axs[1, 6].get_data_ratio())
        if lang == 'en':
            if qty == 'pow':
                axs[1, 6].set_xlabel('Normalized assembly\npower')
            elif qty == 'det':
                axs[1, 6].set_xlabel('Normalized detector\nresponse')
            axs[1, 6].set_ylabel('Probability\ndensity')
        elif lang == 'fr':
            if qty == 'pow':
                axs[1, 6].set_xlabel('Puissance normalisée\nde l\'assemblage')
            elif qty == 'det':
                axs[1, 6].set_xlabel('Réponse normalisée\ndu détecteur')
            axs[1, 6].set_ylabel('Densité\nde\nprobabilité')
        axs[1, 6].tick_params(axis = 'x', labelbottom = True)
        # Add a title and save plot as pdf (vectorized)
        os.system('mkdir -p output_TMC_FSH-BUG_interactions')
        fig.savefig('output_TMC_FSH-BUG_interactions/' + qty + '_' + method
                    + '.pdf',
                    bbox_inches='tight')
        # Clean-up for next plot
        plt.close('all')
#---
#  Plot prior and posterior quantiles of cross sections
#---
EveryXSs = {}
path = '../PyNjoy2016/output/'
directories = glob.glob(path + 'TENDL-2019/*') + glob.glob(path + 'SANDY/*')
# Load any file, just to get the sampled isotopes
lcmfile = loadlcm(links[0])
isotopes = lcmfile['NamIso'].split()
for iso in isotopes:
    print('iso =', iso)
    index = lcmfile[iso]
    # Perform an action only if that isotope has been sampled
    if index != -33:
#        print(iso)
        EveryXSs[iso] = {}
        found = False
        for draglibpath in directories:
            source, isodraglib = draglibpath.rsplit('/', 2)[1:3]
            if isodraglib == iso:
                found = True
                break
        if not found:
            raise Exception('Draglib for ' + iso + ' was not found.')
        print(draglibpath)
        draglib = loadlcm(draglibpath + '/draglib' + iso)
        # Select 2nd temperature (550K), except for H1_H2O (7th, 573.5K)
        if iso == "H1_H2O":
            itemp = 7
        else:
            itemp = 2
        # Initialize the flag indicating presence/absence (resp. True/False) of
        # at least one Autolib data
        AutolibFlag = False
        #---
        #  Establish the list of available random samples, each being contained
        #  in a directory
        #---
        dirnames = [dirnam for dirnam in draglib.keys() if iso in dirnam]
        dirnames.sort()
        #---
        #  Establish the set of available reactions. Threshold reactions may be
        #  present in some random samples, but absent in others (if its cross
        #  sections are too low). The list established here is on a "minimum"
        #  basis, i.e. present in every random samples.
        #---
        reactions = None
        for dirname in dirnames:
            isotopedir = draglib[dirname]
            SubTemp = isotopedir['SUBTMP000' + str(itemp)]
            if not reactions:
                reactions = set(SubTemp.keys())
            else:
                reactions = reactions.intersection(SubTemp.keys())
        #---
        #  Establish the list of interesting reactions among available ones,
        #  for that isotope. The notable reactions set aside are: N3N, N4N,
        #  NNP, N2A, NP, ND
        #---
        wished_reactions = {'NTOT0', 'NG', 'NELAS', 'NFTOT', 'NUSIGF', 'NINEL',
                            'N2N', 'NA'}
        if iso.startswith('U23') or iso.startswith('Zr9'):
            wished_reactions = wished_reactions.difference({'NA'})
        reactions = list(reactions.intersection(wished_reactions))
        reactions.sort()
        labels = {}
        labels['NTOT0'] = "(n,tot)"
        labels['NG'] = r"(n,$\gamma$)"
        labels['NELAS'] = "(n,el)"
        labels['NFTOT'] = "(n,f)"
        labels['NUSIGF'] = r"$\overline{\nu}\times$(n,f)"
        labels['NINEL'] = "(n,inl)"
        labels['N2N'] = "(n,2n)"
        labels['NA'] = r"(n,$\alpha$)"
        labels['NP'] = "(n,p)"
        labels['ND'] = "(n,d)"
        for method in ['BMC', 'BFMC']:
            if method == 'BMC':
                weights = weights_BMC
            elif method == 'BFMC':
                weights = weights_BFMC
            #---
            #  Initialize plot
            #---
            fig, axs = plt.subplots(nrows = len(reactions), sharex = 'all',
                                    figsize = set_size(aspect = 'fullA4'),
                                    gridspec_kw={'hspace': 0.25})
            ireac = -1
            for reaction in reactions:
                ireac = ireac + 1
                print('reaction = ' + reaction)
                # Energy limits (e.g. 99g, 172g, 281g, 295g, 315g, 361g, ...)
                Energies = draglib['ENERGY']
                BaseEnergies = draglib['ENERGY']
                ngrp = np.shape(Energies)[0] - 1
                relstdEnergies = Energies
                # List that shall contain every random samples
                XSs = []
                for dirname in dirnames:
                    irand = int(dirname[-3:])
                    isotopedir = draglib[dirname]
                    # Retrieve the temperature
                    SubTemp = isotopedir['SUBTMP000' + str(itemp)]
                    Temperature = isotopedir['TEMPERATURE']
                    Temperature = str(int(Temperature[itemp - 1]))
                    # Print list of possible reactions
                    if irand == 0:
                        print(SubTemp.keys())
                    # Retrieve cross section
                    XS = SubTemp[reaction]
                    # Combine Autolib data and energy limits (where available) with
                    # XS and energy limits outside of that (where Autolib is
                    # unavailable)
                    if ('BIN-' + reaction) in SubTemp.keys():
                        AutolibFlag = True
                        # Sampling-insensitive data (identical in every random
                        # samples)
                        if irand == 0:
                            AutolibBins = isotopedir['BIN-NFS']
                            BinEnergies = isotopedir['BIN-ENERGY']
                            Autolib_beg = np.min(np.nonzero(AutolibBins))
                            Autolib_end = np.max(np.nonzero(AutolibBins)) + 1
                            Ener_bef, _, Ener_aft = (np.split(Energies,
                                                     [Autolib_beg, Autolib_end + 1]))
                            Energies = np.concatenate((Ener_bef, BinEnergies,
                                                       Ener_aft))
                        XS_bef, _, XS_aft = np.split(XS, [Autolib_beg, Autolib_end])
                        Autolib = SubTemp['BIN-' + reaction]
                        XS = np.concatenate((XS_bef, Autolib, XS_aft))
                    else:
                        # Minor (threshold) reactions may have less than G+1 cross
                        # sections. The missing cross sections shall be considered
                        # as equal to zero. Here, we strip the most thermal energy
                        # limits in order to avoid plotting cross sections equal to
                        # zero.
                        if (len(XS) + 1) < len(Energies):
                            Energies = Energies[:(len(XS) + 1)]
                            relstdEnergies = relstdEnergies[:(len(XS) + 1)]
                    XSs.append(XS)
                    # Deleting progressively to avoid segfaults in lcm module
                    del SubTemp
                    del isotopedir
                nrand = str(len(XSs))
                #---
                #  For some minor (threshold) reactions, we may come up with an
                #  irregular 2D array, because some samples may count more
                #  (non-zero) cross sections than other samples. We shall
                #  homogenize this by reducing to the smallest length.
                #---
                minlength = min(map(len, XSs))
                XSs = [XS[:minlength] for XS in XSs]
                #---
                #  Transform XSs into 2D NumPy array
                #---
                XSs = np.array(XSs)
                #---
                #  Store all the samples of cross sections, for correlations plot
                #---
                nxs = np.shape(XSs)[1]
                if nxs < ngrp:
                    # For threshold reactions, cross sections are zero below that
                    # threshold
                    npad = ((0, 0), (0, ngrp - nxs))
                    EveryXSs[iso][reaction] = np.pad(XSs, pad_width = npad)
                else:
                    EveryXSs[iso][reaction] = XSs
                #---
                #  Add prior quantiles (including median)
                #---
                median = np.quantile(XSs, 0.5, axis=0)
                # Fill-in with light gray between the  5% and 95% quantiles and
                #               dark gray between the 25% and 75% quantiles
                for q1, q2 in [[0.05, 0.95], [0.25, 0.75]]:
                    qXS1 = (np.quantile(XSs, q1, axis=0)/median - 1)*100
                    qXS2 = (np.quantile(XSs, q2, axis=0)/median - 1)*100
                    axs[ireac].fill_between(Energies,
                                            np.append(qXS1, qXS1[-1]),
                                            np.append(qXS2, qXS2[-1]),
                                            linewidth = 0,
                                            step = 'post', color = 'gray',
                                            alpha = 1/(abs(q1 - 0.5)*15 + 1))
                #---
                #  Add posterior weighted quantiles (including median)
                #---
                qn = []
                quantiles = [0.05, 0.25, 0.5, 0.75, 0.95]
                for g in range(0, np.shape(XSs)[1]):
                    qn.append(weighted_quantile(XSs[indices[iso], g],
                                                quantiles,
                                                sample_weight = weights,
                                                values_sorted = False,
                                                old_style = True))
                qn = np.array(qn)
                i = 0
                for q in quantiles:
                    qXS = (qn[:, i]/median - 1)*100
                    axs[ireac].step(Energies, np.append(qXS, qXS[-1]),
                                    where = 'post', linewidth = 0.1,
                                    color = 'r',
                                    alpha = 1/(abs(q - 0.5)*3 + 1))
                    i = i + 1
            #---
            #  Graph glitter
            #---
            ireac = -1
            for reaction in reactions:
                ireac = ireac + 1
                axs[ireac].set_ylabel(labels[reaction])
                # Cross sections are only readable in energy log graphs
                axs[ireac].set_xscale('log')
                # Add energy limits for Autolib domain
                if AutolibFlag:
                    axs[ireac].axvline(BinEnergies[0], linestyle='--',
                                       color='black',
                                       linewidth=0.2)
                    axs[ireac].axvline(BinEnergies[-1], linestyle='--',
                                       color='black',
                                       linewidth=0.2)
                # Restrict the x-axis to data of interest. In particular,
                # remove widely on the left side : the most thermal group is
                # enormous (~3 decades) and of little interest.
                axs[ireac].set_xlim([Energies[-1]/10+Energies[-2]*9/10,
                                     Energies[0]])
                # Use symmetric y limits, capped at 90%
                maxplot = max(np.abs(axs[ireac].get_ylim()))
                maxplot = min(maxplot, 90)*1.05
                maxplot = max(maxplot, 1*1.05) # For O16/B10 unvarying inelastic
                axs[ireac].set_ylim(bottom = -1*maxplot, top = maxplot)
                # Customize y ticks (percent deviation from prior median)
                if len(reactions) <= 3: # Mainly for H in H2O
                    if maxplot < 6:
                        dtick = 2
                    else:
                        dtick = 40
                elif len(reactions) == 4: # Mainly for B10
                    if maxplot < 2:
                        dtick = 1
                    elif maxplot < 4:
                        dtick = 2
                    elif maxplot < 15:
                        dtick = 5
                    else:
                        dtick = 10
                elif len(reactions) == 5: # Mainly for Zr9x
                    if maxplot < 50:
                        dtick = 25
                    else:
                        dtick = 50
                elif len(reactions) == 6: # Mainly for B11, O16, Cr52, Fe5x, Ni58
                    if maxplot < 2:
                        dtick = 1
                    elif maxplot < 10:
                        dtick = 5
                    elif maxplot < 20:
                        dtick = 10
                    elif maxplot < 25:
                        dtick = 15
                    elif maxplot < 50:
                        dtick = 25
                    else:
                        dtick = 50
                elif len(reactions) >= 7: # Mainly for U23x
                    if maxplot < 10:
                        dtick = 5
                    elif maxplot < 20:
                        dtick = 10
                    elif maxplot < 25:
                        dtick = 20
                    elif maxplot < 45:
                        dtick = 25
                    else:
                        dtick = 50
                ticks = [i*dtick for i in range(1, 5)]
                majorticks = list(-np.flip(ticks)) + [0] + ticks
                axs[ireac].yaxis.set_major_locator(ticker.FixedLocator(majorticks))
                ticklabels = [str(i) + '\%' for i in ticks]
                ticklabels = (['-' + i for i in np.flip(ticklabels)] + ['0']
                              + ticklabels)
                axs[ireac].set_yticklabels(ticklabels)
                minorticks = [(i-1/2)*dtick for i in range(1, 5)]
                minorticks = list(-np.flip(minorticks)) + minorticks
                axs[ireac].yaxis.set_minor_locator(ticker.FixedLocator(minorticks))
                # Show all powers of 10 in Energies xaxis, with 10 minor ticks
                locmaj = ticker.LogLocator(base = 10.0, numticks = 24)
                locmin = ticker.LogLocator(base = 10.0, numticks = 24,
                                           subs = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,
                                                   0.8,0.9))
                axs[ireac].xaxis.set_major_locator(locmaj)
                axs[ireac].xaxis.set_minor_locator(locmin)
                axs[ireac].xaxis.set_minor_formatter(ticker.NullFormatter())
                # Add a light grid
                axs[ireac].grid(which='both', alpha=0.3, linewidth=0.1)
            # Set axis label on last subplot
            axs[ireac].set_xlabel('Énergie [eV]')
            # Set xaxis ticks on top instead of bottom
            axs[0].xaxis.set_ticks_position('top')
            # Align labels
            fig.align_ylabels(axs)
            #---
            #  Save plot as pdf (vectorized)
            #---
            os.system('mkdir -p output_TMC_FSH-BUG_interactions')
            fig.savefig('output_TMC_FSH-BUG_interactions/XS_' + iso + '_'
                        + method + '.pdf', bbox_inches='tight')
        del draglib
#---
#  Compute the prior and posterior correlation matrix
#  Ref. : https://tendl.web.psi.ch/bib_rochman/epjp1900495-offprints.pdf
#---
# The most thermal value is not plotted : it takes a lot of space but adds very
# little visual information. Also, switch fast->thermal ordering to
# thermal->fast ordering. Finally, use a log-log scale.
Energies = np.log(np.flip(BaseEnergies)[1:])
ngrp = len(Energies) - 1

# Produce both a global view, with all the isotopes, and also a view zoomed on
# two selected isotopes (uranium 238 and hydrogen)
for view in ['Global', 'Zoom']:
    if view == 'Global':
        isotopes = ['U235', 'U238', 'H1_H2O', 'B10', 'B11', 'O16', 'Cr52',
                    'Fe54', 'Fe56', 'Ni58', 'Zr90', 'Zr91', 'Zr92', 'Zr94',
                    'Zr96']
    elif view == 'Zoom':
        isotopes = ['U238', 'H1_H2O']
    mat = []
    corrlabels = []
    nreaciso = 0
    for iso in isotopes:
        # Matrix correlation is plotted only for some selected (important)
        # reactions
        if iso.startswith('U23'):
            wished_reactions = ['NG', 'NELAS', 'NFTOT', 'NINEL', 'NA']
        elif (iso == 'O16') or (iso == 'B10'):
            wished_reactions = ['NG', 'NELAS', 'NA']
        elif (iso.startswith('Zr9') or (iso == 'Fe54') or (iso == 'Ni58') or
              iso == 'Cr52'):
            wished_reactions = ['NELAS']
        else:
            wished_reactions = ['NG', 'NELAS']
        for reaction in wished_reactions:
            if reaction in EveryXSs[iso].keys():
                mat.append(np.rot90(EveryXSs[iso][reaction][indices[iso], 1:]))
                # Construct correlation labels
                if iso == 'H1_H2O':
                    isolabel = '$^{1}$H'
                else:
                    A = ''.join([s for s in iso if s.isdigit()])
                    elem = ''.join([s for s in iso if s.isalpha()])
                    isolabel = '$^{' + A + '}$' + elem
                corrlabels.append(isolabel + labels[reaction])
                nreaciso = nreaciso + 1
    # Flatten all the dimensions (on isotopes, reactions, energies) except the
    # last one (on random samples)
    mat = np.array(mat)
    mat = mat.reshape(-1, mat.shape[-1])
    #---
    #  Compute and plot this large correlation matrix
    #---
    for method in ['Prior', 'BMC', 'BFMC']:
        # Compute the correlation matrix
        if method == 'Prior':
            Correlation = np.corrcoef(mat)
            prior = Correlation
        else: # Posterior (BMC or BFMC)
            Correlation = weighted_corrcoef(weights, mat)
            # Put back prior's NaNs to posterior
            Correlation = np.where(np.isnan(prior), prior, Correlation)
        # Rotate it, for plotting purposes
        Correlation = np.rot90(Correlation)
        #for delta in [0.5]:
        #for delta in [0.1, 0.2]:
        for delta in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:
            # Initialize the plot
            fig, axs = plt.subplots(nreaciso, nreaciso,
                                    figsize = set_size('correlation'),
                                    gridspec_kw = {'hspace': 0, 'wspace': 0})
            for i in range(nreaciso):
                for j in range(nreaciso):
                    print(i,j)
                    # Select the small piece of correlation matrix that is
                    # plotted at this iteration (one specific reaction of one
                    # specific isotope)
                    subcor = np.rot90(Correlation[i*ngrp:(i+1)*ngrp,
                                                  j*ngrp:(j+1)*ngrp], 3)
                    # This huge amount of energy detail (295g) is reduced, by
                    # grouping very similar values : correlations often have
                    # large homogeneous chunks. This is a very important
                    # saving, which allows this matrix to be actually plotted.
                    subcor, xEnergies, yEnergies = roundcorr(subcor,
                                                             Energies, delta)
                    # Correlation plotting itself
                    pcol = axs[i, j].pcolormesh(xEnergies, yEnergies, subcor,
                                                cmap = 'coolwarm',
                                                linewidth = 0.0,
                                                rasterized = True,
                                                vmin = -1, vmax = 1)
                    # https://stackoverflow.com/questions/27092991/white-lines-in-matplotlibs-pcolor
                    pcol.set_edgecolor('face')
                    # Remove all xticks and yticks...
                    axs[i, j].set_xticks([])
                    axs[i, j].set_yticks([])
                    # ...except for some exceptional cases : example of energy
                    # limits and reaction/isotope labels
                    if (i == 0) and (j == 0) and (nreaciso < 10):
                        axs[i, j].set_xticks([Energies[0], Energies[-1]])
                        axs[i, j].set_xticklabels(['2.5 meV', '19.6 MeV'])
                        axs[i, j].xaxis.tick_top()
                        axs[i, j].set_xlabel('ln($E$)')
                        axs[i, j].xaxis.set_label_position('top')
                    if i == (nreaciso - 1):
                        axs[i, j].set_xlabel(corrlabels[j], rotation = 90)
                    if j == 0:
                        axs[i, j].set_ylabel(corrlabels[nreaciso-1-i],
                                             rotation = 0,
                                             va = 'center',
                                             labelpad = 30.0)
            # Add a colorbar below. Ref. :
            # https://stackoverflow.com/questions/13784201/how-to-have-one-colorbar-for-all-subplots
            fig.subplots_adjust(bottom=0.05)
            # [0.15 ; -0.12] : Coordinates of lower-left corner of the
            #                  colorbar. Here, these coordinates indicate the
            #                  lower-left corner of the figure.
            # 0.7  : Relative width of the colorbar
            # 0.02 : Relative height of the colobar
            cbar_ax = fig.add_axes([0.15, -0.12, 0.7, 0.02])
            # Colorbar must have discrete colors, in order to faithfully
            # represent the plotted data (rounded values)
            halfscale = np.arange(delta/2, 1, delta)
            halfscale = np.concatenate((halfscale, np.ones(1)))
            scale = np.concatenate((-np.flip(halfscale), halfscale))
            cbar = fig.colorbar(pcol, cax = cbar_ax, orientation = 'horizontal',
                                boundaries = scale, spacing = 'proportional',
                                drawedges = True)
            # Set thinner limits between each blocks of homogeneous discrete
            # color
            cbar.dividers.set_linewidth(0.25)
            # Set ticks independently of colors
            if math.isclose(delta, 0.5) or math.isclose(delta, 0.3):
                ticks = np.concatenate((-np.flip(halfscale), [0], halfscale))
                cbar.set_ticks(ticks)
                ticklabels = halfscale[:-1]
                ticklabels = (['-1']
                              + ["-"
                              + "{:.2f}".format(k) for k in np.flip(ticklabels)]
                              + ['0']
                              + ["{:.2f}".format(k) for k in ticklabels] + ['1'])
            else:
                # Tick every 0.2
                cbar.set_ticks(np.linspace(-1, 1, num = 11))
                ticklabels = np.linspace(0, 1, num = 6)[1:-1]
                ticklabels = (['-1']
                              + ["-"
                              + "{:.1f}".format(k) for k in np.flip(ticklabels)]
                              + ['0']
                              + ["{:.1f}".format(k) for k in ticklabels] + ['1'])
            cbar.ax.set_xticklabels(ticklabels)
            # Change label size in the colorbar
            cbar.ax.tick_params(labelsize = 12)
            # Save plot as pdf (vectorized, except for the color dots)
            fig.savefig('output_TMC_FSH-BUG_interactions/Correlation' + view
                        + '_' + method + '_' + str(int(delta*100)) + '.pdf',
                        bbox_inches='tight', dpi = 600)
            # Clean-up for next plot
            plt.close('all')
            del fig
            del axs
#---
#  Convergence plot: standard deviation of detector reading in one assembly (H1)
#---
# Retrieve detector reading in assembly localized in H1
stat = responsesfTMC[:, 0, 0]
# Create a new plot
fig, ax = plt.subplots(figsize = set_size())
for method in ['Prior', 'BFMC', 'BMC']:
    # Prepare legend labels and retrieve weights
    if method == 'BMC':
        label = 'A posteriori (BMC)'
        weights = weights_BMC
    elif method == 'BFMC':
        label = 'A posteriori (BFMC)'
        weights = weights_BFMC
    else:
        label = 'A priori'
        weights = np.ones_like(stat)/len(stat)
    # Compute cumulative detector reading
    nrand = len(stat)
    xaxis = np.arange(0, nrand)
    yaxis = np.zeros_like(stat)
    for i in xaxis:
        yaxis[i] = weighted_relstd(stat[:i + 1], weights[:i + 1])
    # Add data
    ax.plot(xaxis, yaxis, label = label)
# Graph glitter
plt.legend(framealpha = 1.0, loc = 'lower right')
ax.yaxis.set_major_formatter(ticker.PercentFormatter())
ax.grid()
ax.set_ylim(bottom = 0)
ax.set_xlim(left = 0, right = nrand)
ax.set_xlabel('Nombre d\'échantillons')
ax.set_ylabel('Écart-type relatif de la réponse\ndu détecteur localisé en H1')
# Save plot as pdf (vectorized)
fig.savefig('output_TMC_FSH-BUG_interactions/Conv_det.pdf',
            bbox_inches='tight')
# Clean-up for next plot
plt.close('all')
#---
#  Convergence plot: 5 quantiles of 238U(n,inl) cross section at 3 MeV
#---
# Retrieve inelastic cross section of U238 at 3 MeV
iso = 'U238'
reac = 'NINEL'
# Filter only highest energies
stat = EveryXSs[iso][reac][:, BaseEnergies[:-1] > 3e+6]
# Retrieve the 2nd lowest of these energies, for nrand samples
stat = stat[indices[iso], -2]
# Create a new plot
fig, ax = plt.subplots(figsize = set_size())
for method in ['Prior', 'BFMC', 'BMC']:
    # Prepare legend labels and retrieve weights
    if method == 'BMC':
        label = 'A posteriori (BMC)'
        color = '#2ca02c' # Green
        weights = weights_BMC
    elif method == 'BFMC':
        label = 'A posteriori (BFMC)'
        color = '#ff7f0e' # Orange
        weights = weights_BFMC
    else:
        label = 'A priori'
        color = '#1f77b4' # Blue
        weights = np.ones_like(stat)/len(stat)
    # Compute cumulative quantiles
    quantiles = [0.05, 0.25, 0.5, 0.75, 0.95]
    nrand = len(stat)
    xaxis = np.arange(0, nrand)
    yaxis = np.zeros((len(stat), len(quantiles)))
    for i in xaxis:
        yaxis[i, :] = weighted_quantile(stat[:i + 1],
                                        quantiles,
                                        sample_weight = weights[:i + 1],
                                        values_sorted = False,
                                        old_style = True)
    for q in range(0, len(quantiles)):
        # Add data
        if q == 2:
            ax.plot(xaxis, yaxis[:, q], color = color, label = label)
        else:
            if (q == 1) or (q == 3):
                linestyle = '--'
            else:
                linestyle = ':'
            ax.plot(xaxis, yaxis[:, q], color = color, linestyle = linestyle,
                    alpha = 1/(abs(quantiles[q] - 0.5)*3 + 1))
# Graph glitter
plt.legend(framealpha = 1.0, loc = 'lower right')
ax.grid()
ax.set_xlim(left = 0, right = nrand)
ax.set_xlabel('Nombre d\'échantillons')
ax.set_ylabel('Section efficace de diffusion inélastique'
              + '\nde l\'uranium 238 à 3\,MeV [b]')
# Save plot as pdf (vectorized)
fig.savefig('output_TMC_FSH-BUG_interactions/Conv_U8inl.pdf',
            bbox_inches='tight')
# Clean-up for next plot
plt.close('all')
#---
#  Convergence plot: correlation betwen
#  * 238U(n,inl) cross section at 3 MeV, and
#  * H1(n,n) cross section at 7 MeV
#---
# Retrieve elastic cross section of hydrogen at 7 MeV
iso = 'H1_H2O'
reac = 'NELAS'
# Filter only highest energies
stat2 = EveryXSs[iso][reac][:, BaseEnergies[:-1] > 7e+6]
# Retrieve the 2nd lowest of these energies, for nrand samples
stat2 = stat2[indices[iso], -2]
# Create a new plot
fig, ax = plt.subplots(figsize = set_size())
for method in ['Prior', 'BFMC', 'BMC']:
    # Prepare legend labels and retrieve weights
    if method == 'BMC':
        label = 'A posteriori (BMC)'
        weights = weights_BMC
    elif method == 'BFMC':
        label = 'A posteriori (BFMC)'
        weights = weights_BFMC
    else:
        label = 'A priori'
        weights = np.ones_like(stat)/len(stat)
    # Compute cumulative detector reading
    nrand = len(stat)
    xaxis = np.arange(1, nrand)
    Correl = np.zeros_like(stat)
    for i in xaxis:
        Correl[i] = weighted_corrcoef(weights[:i + 1],
                                      stat[:i + 1],
                                      stat2[:i + 1])[0, 1]
    # Add data
    ax.plot(xaxis, Correl[1:], label = label)
# Graph glitter
plt.legend(framealpha = 1.0, loc = 'lower right')
ax.grid()
ax.set_xlim(left = 0, right = nrand)
ax.set_ylim(bottom = -1)
ax.set_xlabel('Nombre d\'échantillons')
ax.set_ylabel('Corrélation')
# Save plot as pdf (vectorized)
fig.savefig('output_TMC_FSH-BUG_interactions/Conv_correl_H1_U238.pdf',
            bbox_inches='tight')
# Clean-up for next plot
plt.close('all')
