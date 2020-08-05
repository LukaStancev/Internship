#
#  Plotting LCM file contents
#  Usage : python3 Plot.py
#

# Imports
import lcm
import numpy as np
import os

# Retrieve data
os.system("ln -s ../Drakkar/Linux_x86_64/_PowerD.ascii .")
my_lcm=lcm.new('LCM_INP','PowerD.ascii')
os.system("rm -f _PowerD.ascii")

# Data processing
my_lcm._impx=3
power=my_lcm['POWER-CHAN'] # Retrieve
mean=sum(power)/len(power) # Compute mean power
power=power/mean # Normalize
print(power)

#import matplotlib.pyplot as plt
##fig, ax = plt.subplots()
##ax.plot([1, 2, 3, 4], [1, 4, 2, 3])
#
#plt.plot([1, 2, 3, 4], [1, 4, 2, 3])
#plt.savefig('myfig')

print("Plotting completed")
