#  Execute D2S on every available assembly geometries
#  Usage   : python3 Assemblies_kinf.py
#  Date    : 1/2021
import glob
import os
from D2S import InfiniteLattice

for geofilepath in glob.glob('../geo_compo/Almaraz/UOX*.geo'):
    filepath = geofilepath[:-4]
    print(filepath)
    InfiniteLattice(filepath)

os.system('mv *.sss2 ../Assemblies_kinf/.')
