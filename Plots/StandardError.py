# Plot standard errors (SE) for standard deviation (SD), skewness (SK) and
# excess kurtosis (KU) for n between 10 and 10000
# Author : V. Salino, 06/2021

from BasicFunctions import *
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import copy
plt.rcParams.update(tex_fonts())
# Managing the language
lang = 'fr' # fr/en
if lang == 'en':
    xlabel = 'Number of independent samples'
    ylabel1 = ('Standard error of standard deviation (unitless),\n'
               + 'for a normal distribution')
    ylabel2 = ('Standard errors of skewness and excess kurtosis (unitless),\n'
               + 'for any distribution')
    title = ('Standard errors on standard deviation (blue, left axis),\n'
             + 'skewness and excess kurtosis (green, right axis)')
    legend = ['Standard deviation','Skewness', 'Excess kurtosis']
elif lang == 'fr':
    xlabel = 'Nombre d\'échantillons indépendants'
    ylabel1 = ('Erreur-type (sans unité) de l\'écart-type,\n'
               + 'pour une distribution normale')
    ylabel2 = ('Erreurs-types (sans unité) de l\'asymétrie et l\'excès de '
               +'kurtosis,\npour n\'importe quelle distribution de '
               + 'probabilité')
    title = ''
    legend = ['Écart-type', 'Asymétrie', 'Excès de kurtosis']
# Initialize the figure and the number of points
fig, ax1 = plt.subplots(figsize = set_size())
n = np.arange(10,10001)
#---
#   Standard deviation (left axis)
#---
line1 = ax1.plot(n, SE_SD_Normal(n)*100)
ax1.set(xlabel = xlabel, ylabel = ylabel1)
if title:
    ax1.set(title = title)
# That plot is much easier to read on a log-scale
ax1.set_xscale('log')
# Restrict to the interesting part
ax1.set_xlim([10, 10000])
ax1.set_ylim([0, 24])
# Add mid-decade xticks and scalar formatting
ax1.set_xticks([10, 50, 100, 300, 1000, 10000])
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# Full grid, to allow reading values on the plot
ax1.grid(which = 'major', linewidth = 1.3)
ax1.grid(which = 'minor', linewidth = 0.8)
# Major tick every 2% and minor tick every %
ax1.yaxis.set_major_locator(MultipleLocator(2))
ax1.yaxis.set_minor_locator(MultipleLocator(1))
# Add percent sign after yticks
labels = ax1.get_yticks()
percentlabels = (['invisible', 0]
                 + ['{:1.0f}'.format(i) + '\%' for i in labels[2:]])
ax1.set_yticklabels(percentlabels)
# Dealing with doubled axes
ax1.tick_params(axis='y', color='C0', labelcolor='C0')
ax1.yaxis.label.set_color('C0')
ax1.spines['left']. set_visible(False)
ax1.spines['right']. set_visible(False)
#---
#  Skewness and kurtosis (right axis)
#---
ax2 = ax1.twinx()
line2 = ax2.plot(n, SE_SK(n), color = 'C2')
line3 = ax2.plot(n, SE_KU(n), color = 'C2', linestyle="--")
ax2.set_ylabel(ylabel2)
# Restrict to the interesting part
ax2.set_ylim([0, 1.2])
# Major tick every 0.1
ax2.yaxis.set_major_locator(MultipleLocator(0.1))
# Dealing with doubled axes
ax2.tick_params(axis='y', color='C2', labelcolor='C2')
ax2.yaxis.label.set_color('C2')
ax2.spines['right'].set_color('C2')
ax2.spines['left'].set_color('C0')
# Add a legend
lines = line1 + line2 + line3
leg = ax2.legend(lines, legend, handlelength = 2.5, framealpha = 1)
# Set the linewidth of each legend object
for legobj in leg.legendHandles:
        legobj.set_linewidth(1.8)
#---
#  Save plot as pdf (vectorized)
#---
fig.savefig('StandardError.pdf', bbox_inches = 'tight')
