#  Execute D2S on every full core geometries
#  Usage   : python3 Full_core.py
#  Date    : 1/2021
import glob
import os
from D2S import FullCore

# Full core geometry, suited for every possibility of control rod insertions
filepath = '../geo_compo/FullCore.geo'

# All rods out
Assemblies = {}
Assemblies[4] = '195_None'
Assemblies[5] = '255_Py8'
Assemblies[6] = '255_Py12'
Assemblies[7] = '310_None'
Assemblies[8] = '310_Py12'
Assemblies[9] = '195_None'
Assemblies[10] = '195_None'
CB = 1206 # ppm
FullCore(filepath, Assemblies, CB)

# D control rod bank inserted
Assemblies[10] = '195_AIC'
CB = 1084 # ppm
FullCore(filepath, Assemblies, CB)

# C and D control rod banks inserted
Assemblies[9] = '195_AIC'
CB = 960 # ppm
FullCore(filepath, Assemblies, CB)

os.system('mv *.sss2 ../.')
