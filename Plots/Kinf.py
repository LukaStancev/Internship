#
#  Plotting K-infinite difference between Dragon and Serpent
#  Usage : python3 Kinf.py
#  Author : V. Salino (IRSN), 02/2021
#

# Imports
# faulthandler for segmentation faults. Left here due to a heisenbug in lcm
import faulthandler; faulthandler.enable()
import lcm
import serpentTools
from serpentTools.settings import rc
import numpy as np
from scipy import interpolate
import os
import glob
import re
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#---
#  Retrieve Serpent k-inf
#---
# Dictionary intended to receive Serpent data (plus/minus standard deviation)
Serpent = {}
Serpentpstd = {}
Serpentmstd = {}
# Loop for all Serpent
respaths = glob.glob('../Serpent/Assemblies_kinf/UOX*.sss2_res.m')
respaths.sort()
for respath in respaths:
    # Retrieve Serpent k-inf absorption-based estimator and its standard
    # deviation
    res = serpentTools.read(respath)
    kinf = res.resdata['absKinf']
    # Retrieve type (enrichment, control rods insertion) and boron
    # concentration
    filename = respath.split('/')[-1].split('.')[0]
    Enrich, Rods, CB = filename.split('_')
    CB = int(re.sub('[^0-9]', '', CB))
    Type = re.sub('[^0-9]', '', Enrich) + "_" + Rods
    # Store all this data
    if Type not in Serpent:
        Serpent[Type] = []
        Serpentpstd[Type] = []
        Serpentmstd[Type] = []
    Serpent[Type].append([CB, kinf[0]])
    Serpentpstd[Type].append([CB, kinf[0]+kinf[1]])
    Serpentmstd[Type].append([CB, kinf[0]-kinf[1]])
# For each type, get into a 2D numpy array and sort by boron
for Type in Serpent:
    Serpent[Type] = np.array(Serpent[Type])
    Serpent[Type] = Serpent[Type][Serpent[Type][:,0].argsort()]
    Serpentpstd[Type] = np.array(Serpentpstd[Type])
    Serpentmstd[Type] = np.array(Serpentmstd[Type])
    Serpentpstd[Type] = Serpentpstd[Type][Serpentpstd[Type][:,0].argsort()]
    Serpentmstd[Type] = Serpentmstd[Type][Serpentmstd[Type][:,0].argsort()]
#---
#  Retrieve Dragon k-inf
#---
# Dictionary intended to receive Dragon data
Dragon = {}
for Multicompofile in glob.glob('../Drakkar/Linux_x86_64/_UOX*.ascii'):
    Type = Multicompofile.split('/')[-1][4:].split('.')[0]
    os.system('ln -s ' + Multicompofile + ' _Multicompo.ascii')
    Multicompo = lcm.new('LCM_INP', 'Multicompo.ascii')
    os.system('rm -f _Multicompo.ascii')
    #---
    #  See 'A Description of the DRAGON and TRIVAC Version5 Data Structures'
    #  (IGE351), Contents of a /saphyb/ directory
    #---
    # Retrieve parameter names
    PARKEY = Multicompo['paramdescrip']['PARKEY']
    PARKEY = [PARKEY[i:i + 4] for i in range(0, len(PARKEY), 4)]
    # Retrieve boron concentrations and control rods insertions
    for pval in sorted(Multicompo['paramvaleurs'].keys()):
        ipar = int(pval.split()[1])
        if PARKEY[ipar-1] == 'CBOR':
            CBOR = Multicompo['paramvaleurs'][pval]
            iparcb = ipar
        elif PARKEY[ipar-1] == 'BARR':
            BARR = Multicompo['paramvaleurs'][pval]
            iparcr = ipar
        else:
            if len(Multicompo['paramvaleurs'][pval]) != 1:
                raise Exception('This script has been designed to take as '
                                + 'input a library with boron and control rods'
                                + ' variations only. However, '
                                + PARKEY[ipar-1] + ' has abnormally more than '
                                + 'one value : '
                                + str(Multicompo['paramvaleurs'][pval]))
    # Retrieve tree parameters
    NVP = Multicompo['DIMSAP'][16]
    NCALS = Multicompo['DIMSAP'][18]
    DEBARB = Multicompo['paramarbre']['DEBARB']
    ARBVAL = Multicompo['paramarbre']['ARBVAL']
    npar = int(Multicompo['paramdescrip']['NPAR'])
    # Initialize data to be plotted
    kinf = np.zeros(NCALS)
    CB = np.zeros(NCALS)
    CR = np.zeros(NCALS)
    # Go through the elementary calculations
    for calc in [Dir for Dir in Multicompo.keys() if Dir.startswith('calc')]:
        # Retrieve parameter m-uplet
        ical = int(calc.split()[1])
        for i in range(NVP - NCALS + 1, NVP + 1):
            if DEBARB[i] == ical:
                i0 = i
        muplet = np.zeros(npar)
        muplet[npar - 1] = ARBVAL[i0 - 1]
        for ipar in range(npar-1, 0, -1):
            i0 = np.where(DEBARB[:NVP - NCALS + 1] > i0)[0][0]
            muplet[ipar - 1] = ARBVAL[i0 - 1]
        # Retrieve data to be plotted
        kinf[ical-1] = Multicompo[calc]['divers']['VALDIV'][1]
        CB[ical-1] = CBOR[int(muplet[iparcb-1] - 1)]
        CR[ical-1] = BARR[int(muplet[iparcr-1] - 1)]
    # Store this data in a dictionnary
    if '_Py' not in Type:
        if (Type + '_None') not in Dragon:
            Dragon[Type + '_None'] = []
        if (Type + '_AIC') not in Dragon:
            Dragon[Type + '_AIC'] = []
        indices = [i for i,x in enumerate(CR == 1) if x]
        Dragon[Type + '_None'].append(CB[indices])
        Dragon[Type + '_None'].append(kinf[indices])
        indices = [i for i,x in enumerate(CR == 2) if x]
        Dragon[Type + '_AIC'].append(CB[indices])
        Dragon[Type + '_AIC'].append(kinf[indices])
    else:
        if Type not in Dragon:
            Dragon[Type] = []
        Dragon[Type].append(CB)
        Dragon[Type].append(kinf)
# For each type, get into a 2D numpy array and sort by boron
for Type in Dragon:
    Dragon[Type] = np.array(Dragon[Type])
    Dragon[Type] = Dragon[Type].transpose()
#---
#  Compute k-inf discrepancy (with natural logarithm)
#---
Discrepancies = {}
Discrepanciespstd = {}
Discrepanciesmstd = {}
for Type in Serpent:
    print(Type)
    f = interpolate.interp1d(Dragon[Type][:,0], Dragon[Type][:,1],
                             kind = 'quadratic')
    DragonInterp = f(Serpent[Type][:,0])
    Discrepancies[Type] = np.log(Serpent[Type][:,1]/DragonInterp)*1E5
    Discrepanciespstd[Type] = np.log(Serpentpstd[Type][:,1]/DragonInterp)*1E5
    Discrepanciesmstd[Type] = np.log(Serpentmstd[Type][:,1]/DragonInterp)*1E5
    print(Discrepanciesmstd[Type])
    print(Discrepancies[Type])
    print(Discrepanciespstd[Type])
#---
#  Produce a bar plot
#---
nbBor = len(list(Discrepancies.values())[0])
w = 0.75
dimw = w/nbBor
fig, ax = plt.subplots()
x = np.arange(len(Discrepancies))
for i in range(nbBor):
    y = []
    ypstd = []
    ymstd = []
    Types = []
    xlabels = []
    SortedDiscr = sorted(Discrepancies.items(), key=lambda x: x[1][1])
    for Type, Discrepancy in SortedDiscr:
        # Get mean discrepancy
        y.append(Discrepancy[i])
        # Construct legend from 'Type'. For example: '195_None', '255_Py8', ...
        Enrichment, Rods = Type.split('_')
        Enrichment = Enrichment[0] + '.' + Enrichment[1:] + '%'
        if Rods == 'None':
            Rods = ''
        elif Rods == 'AIC':
            Rods = '20 AIC rods'
        elif Rods.startswith('Py'):
            Rods = Rods[2:] + ' Pyrex\nrods'
        xlabels.append(Enrichment + '\n' + Rods)
        # Get discrepancy for +/- 1 sigma
        ypstd.append(Discrepanciespstd[Type][i])
        ymstd.append(Discrepanciesmstd[Type][i])
    ax.bar(x + i * dimw - 0.25, y, dimw,
           yerr = (np.array(y) - np.array(ymstd),
                   np.array(ypstd) - np.array(y)),
           error_kw=dict(capsize = w*10, capthick = w*2/3, lw=w*2/3),
           label = str(int(Serpent[Type][i,0])) + ' ppm of boron')
#---
#  Graph glitter
#---
plt.xticks(x, xlabels)
ax.set_xlabel('Assembly ' + r'$^{235}\mathrm{U}$'
              + ' enrichment and its eventual absorbing rods')
ax.set_ylabel(r'$\mathrm{ln}\left(\frac{'
              r'k^{\mathrm{Serpent}}_{\mathrm{infinite}}}{'
              r'k^{\mathrm{Dragon}}_{\mathrm{infinite}}}'
              r'\right)\times10^5$', fontsize = 14)
ax.set_title('Discrepancies of infinite multiplication factors\n'
             + 'between Serpent and Dragon, in pcm,\n'
             + 'for the assemblies in Tihange at the first startup.\n'
             + 'Error bars correspond to the Monte-Carlo $1\sigma$ '
             + 'uncertainty.')
plt.legend()
# Add minor ticks
ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))
# Tight layout prevents a label from exceeding the figure frame
fig.tight_layout()
fig.savefig('Kinf.pdf')
print("Plotting completed")
