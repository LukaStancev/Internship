#
#  Compute Pearson correlation coefficient (squared) between:
#  * power of the central assembly, and
#  * hydrogen (H1) cross sections
#  Usage : python3 Pearson.py
#  Author : V. Salino (IRSN), 02/2021
#

# Imports
import lcm
import numpy as np
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import pearsonr
import math
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from BasicFunctions import *
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en

# Retrieve H1 cross sections, through (very basic) ENDF parsing
from LoadH1H2O import *
os.system('./ParseH1H2O.sh')
energies, XSs = LoadH1H2O()
# Initialize a gaussian vector (used for legends)
gaussian = np.random.randn(300)
# Declare control rod insertions and its full names
controlrods = ['CD', 'D', 'ARO']
text = {}
text['ARO'] = 'all rods out'
text['D'] = 'D rod bank inserted'
text['CD'] = 'C and D rod banks inserted'
# Create a dictionnary to store Pearson Correlation Coefficient and its pvalue
# between the central assembly and the H1_H2O cross section
PCC = {}
pvalue = {}
for controlrod in controlrods:
    print('controlrod = ' + controlrod)
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
    isamples = []
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
            if (int(ResultFile[iso]) != -33) and (iso == 'H1_H2O'):
                powers[iso].append(power)
                # Also recover the sample number, to order them and compute
                # correlations
                isamples.append(int(ResultFile['H1_H2O']))
        del ResultFile
    # Update 'isotopes' list so that only those isotopes that were indeed
    # randomly sampled can be the subject of loops to come
    isotopes_used = []
    for iso in isotopes:
        if len(powers[iso]) != 0:
            isotopes_used.append(iso)
    isotopes = isotopes_used
    # Order them, to facilitate the calculation of correlations
    for iso in isotopes:
        print('current :' + iso)
        tmp = []
        for (a,b) in sorted(zip(isamples, powers[iso])):
            tmp.append(b)
        powers[iso] = tmp
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
    #  Compute Pearson's correlation coefficient (PCC) for the central assembly
    #---
    PCC[controlrod] = {}
    pvalue[controlrod] = {}
    for reaction in XSs:
        print(reaction)
        PCC[controlrod][reaction] = np.zeros_like(energies[reaction])
        pvalue[controlrod][reaction] = np.zeros_like(energies[reaction])
        for i in np.arange(0, len(energies[reaction])):
            # Focus on the central assembly
            PCC[controlrod][reaction][i], pvalue[controlrod][reaction][i] = pearsonr(XSs[reaction][i, :], Powers2D['H1_H2O'][:, 7, 0])
#---
#  Correlations plots
#---
for controlrod in ['CD', 'D', 'ARO']:
    fig, ax = plt.subplots()
    for reaction in XSs:
        if reaction == 'nelastic':
            txt = '(n,el)'
        elif reaction == 'ngamma':
            txt = r'(n,$\gamma$)'
        else:
            txt = 'void'
        ax.step(energies[reaction][:-1], PCC[controlrod][reaction][:-1]**2,
                where = 'post', label = txt)
    #---
     #  Graph glitter
    #---
    if lang == 'en':
      ax.set_title('Correlations between the uncertainty of the central '
                   + 'assembly\npower and the elastic and capture cross '
                   + ' sections of hydrogen,\n' + text[controlrod])
      ax.set_xlabel('Energy [eV]')
      ax.set_ylabel('Pearson correlation coefficient, squared')
    elif lang == 'fr':
      ax.set_xlabel('Énergie [eV]')
      ax.set_ylabel('Coefficient de corrélation de Pearson au carré')
    ax.set_xscale('log')
    ax.legend(loc = 'center right')
    ax.set_xlim(left = min(energies[reaction]),
                right = max(energies[reaction][:-1]))
    ax.set_ylim(bottom = 0, top = 1)
    ax.grid(which='both', alpha=0.2, linewidth=0.1)
    #---
    #  Save plot as pdf (vectorized)
    #---
    os.system('mkdir -p output_Pearson')
    fig.savefig('output_Pearson/' + controlrod + '.pdf', bbox_inches='tight')
    #---
    #  Clean-up
    #---
    plt.close('all')

