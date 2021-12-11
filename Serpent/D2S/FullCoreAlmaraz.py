#  Execute D2S on Almaraz-2 full core geometry
#  Usage   : python3 Full_core.py
#  Date    : 1/2021
import glob
import os
from D2S import FullCore

# Full core geometry, suited for every possibility of control rod insertions
filepath = '../geo_compo/Almaraz/FullCore.geo'

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
CB = 1280 # ppm
FullCore(filepath, Assemblies, CB)

os.system('mv *.sss2 ../FullCore/.')
