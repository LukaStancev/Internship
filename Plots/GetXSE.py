# Imports
# faulthandler for segmentation faults. Left here due to a heisenbug in lcm
import faulthandler; faulthandler.enable()
import lcm
import numpy as np
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from BasicFunctions import *

# Function that lists all available isotope datas
def GetISO():
	 #graphtype = 'samples'
 path = '../PyNjoy2016/output/'
 directories = glob.glob(path + 'SANDY/*')
 minsig = 0.01
 #list initialization
 iso_list = []  
 for draglibpath in directories:
         # Retrieve isotope name as given in directories' names and its source
         source, iso = draglibpath.rsplit('/', 2)[1:3]
         # adding elements to the list
         iso_list.append(iso)
 for iso in iso_list:
         print('isotopes = ', iso)
 return 

# Function that reads Draglib data for certain choosen isotope
def GetXSE(choice):

 #graphtype = 'samples'
 path = '../PyNjoy2016/output/'
 directories = glob.glob(path + 'SANDY/*')
 minsig = 0.01
 #list initialization
 iso_list = []  
 for draglibpath in directories:
         # Retrieve isotope name as given in directories' names and its source
         source, iso = draglibpath.rsplit('/', 2)[1:3]
         # adding elements to the list
         iso_list.append(iso)
         
 if choice in iso_list:
       iso = str(choice)
       draglibpath = draglibpath.replace(iso_list[-1], iso)
       print("You selected", iso)
 else:
       print("Not a good value")
#---
#  Create link toward Draglib file and load data (most time spent here)
#---
 print('ln -s ' + draglibpath + '/draglib' + iso + ' _draglib' + iso
               + '.txt')
 command = ('ln -s ' + draglibpath + '/draglib' + iso + ' _draglib' + iso
               + '.txt')
 os.system(command)
 draglib = lcm.new('LCM_INP', 'draglib' + iso + '.txt')
 os.system('rm -f _draglib' + iso + '.txt')
 draglib._impx = 0 

# Select 2nd temperature (550K), except for H1_H2O (7th, 573.5K)
 if iso == "H1_H2O":
   itemp = 7
 else:
   itemp = 2

 # Initialize the flag indicating presence/absence (resp. True/False) of at
 # least one Autolib data
 AutolibFlag = False
 #---
 #  Establish the list of available random samples, each being contained in
 #  a directory
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
#  Establish the list of interesting reactions among available ones, for
#  that isotope. The notable reactions set aside are: NG, NELAS, NFTOT,
#  N3N, N4N, NNP, N2A
#---
 wished_reactions = {'NTOT0', 'NUSIGF', 'NINEL', 'N2N', 'NA', 'NP', 'ND'}
 if iso.startswith('U23'):
    wished_reactions = wished_reactions.difference({'NA', 'NP', 'ND'})
 elif iso.startswith('Zr9'):
    wished_reactions = wished_reactions.difference({'NA'})
 reactions = list(reactions.intersection(wished_reactions))
 reactions.sort()
 labels = {}
 labels['NTOT0'] = "(n,tot)"
 labels['NUSIGF'] = r"$\overline{\nu}\times$(n,f)"
 labels['NINEL'] = "(n,inl)"
 labels['N2N'] = "(n,2n)"
 labels['NA'] = r"(n,$\alpha$)"
 labels['NP'] = "(n,p)"
 labels['ND'] = "(n,d)"
 for reaction in reactions:
        print('reaction = ' + reaction)
        
#---
#  Retrieve data that shall be plotted
#---
# Energy limits (e.g. 99g, 172g, 281g, 295g, 315g, 361g, ...)
 Energies = draglib['ENERGY']
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
            # Combine Autolib data and energy limits (where available) with XS
            # and energy limits outside of that (where Autolib is unavailable)
            if ('BIN-' + reaction) in SubTemp.keys():
                AutolibFlag = True
                # Sampling-insensitive data (identical in every random samples)
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
                # sections. The missing cross sections shall be considered as
                # equal to zero. Here, we strip the most thermal energy limits
                # in order to avoid plotting cross sections equal to zero.
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
#  irregular 2D array, because some samples may count more (non-zero)
#  cross sections than other samples. We shall homogenize this by
#  reducing to the smallest length.
#---
 minlength = min(map(len, XSs))
 XSs = [XS[:minlength] for XS in XSs]
#---
#  Transform XSs into 2D NumPy arrays and then compute mean and
#  relative standard deviation (in %) over every samples
#---
 XSs = np.array(XSs)
 XS = np.mean(XSs, axis = 0)
 relstdXS = np.std(XSs, axis = 0, ddof = 1)/XS*100
#---
#  Keep relative standard deviation only when the cross section is
#  larger than minsig b (for threshold reactions)
#---
 if XS[0] > minsig:
            largeXS = np.where(XS > minsig)[0]
            limit = len(np.intersect1d(largeXS, np.arange(0, len(largeXS))))
            relstdXS = relstdXS[:limit]
            relstdEnergies = Energies[:(len(relstdXS) + 1)]
            np.savetxt('XSGETXSE.txt', XSs)
 del draglib
 # Returning the results
 return relstdEnergies, relstdXS, Energies, XS, reactions, iso, iso_list 
