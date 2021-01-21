#  Execute D2S on every available geometries
#  Usage   : python3 Tihange.py
#  Date    : 1/2021
import glob
import os
from D2S import D2S

for geofilepath in glob.glob('../geo_compo/UOX*.geo'):
    filepath = geofilepath[:-4]
    print(filepath)
    D2S(filepath)

os.system('mv *.sss2 ../.')
