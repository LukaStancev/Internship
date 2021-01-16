#  Execute D2S on every available geometries
#  Usage   : python3 Tihange.py
#  Date    : 1/2021
import glob
import os
from D2S import D2S

for geofilename in glob.glob('_UOX*.geo'):
    print(geofilename)
    D2S(geofilename)

os.system('mv *.sss2 ../.')
