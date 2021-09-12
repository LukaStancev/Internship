#  Execute D2S on every full core geometries of Tihange-1
#  Usage   : python3 Full_core.py
#  Date    : 1/2021
import glob
import os
from D2S import FullCore

# Full core geometry, suited for every possibility of control rod insertions
filepath = '../geo_compo/Fessenheim-Bugey/FullCore.geo'

# All rods out
Assemblies = {}
Assemblies[4] = '210_None'
Assemblies[5] = '310_None'
Assemblies[6] = '260_Py12a'
Assemblies[7] = '260_Py16a'
Assemblies[8] = '260_Py20a'
Assemblies[9] = '310_Py12a'
Assemblies[10] = '310_Py16a'
Assemblies[11] = '210_AIC'
CB = 1325 # ppm
FullCore(filepath, Assemblies, CB)

os.system('mv *.sss2 ../FullCore/.')
