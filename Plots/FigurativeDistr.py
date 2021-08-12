from BasicFunctions import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('pdf')
import numpy as np
import mpmath
import math
plt.rcParams.update(tex_fonts())
lang = 'fr' # fr/en

# Gaussian function
def gaussian(x, mu, sigma):
    y = np.zeros_like(x)
    for i in range(0, len(x)):
        y[i] = ( math.exp(-(((x[i] - mu)/sigma)**2)/2) /
                (sigma * math.sqrt(2*math.pi)) )
    return y
# Hyperbolic secant distribution
# https://fr.wikipedia.org/wiki/Loi_s%C3%A9cante_hyperbolique
def sech(x):
    y = np.zeros_like(x)
    for i in range(0, len(x)):
        y[i] = mpmath.sech(math.pi/2*x[i])/2
    return y
# Wigner semicircle distribution
# https://fr.wikipedia.org/wiki/Loi_du_demi-cercle
def wigner(x, sigma):
    y = np.zeros_like(x)
    R = 2*sigma
    for i in range(0, len(x)):
        if -R < x[i] < R:
            y[i] = 2/(math.pi*R**2)*math.sqrt(R**2 - x[i]**2)
    return y

# Focus on right tail
x = np.linspace(1.9, 4, 1000)
# Same mean and standard deviation for all the distributions
mu = 0
sigma = 1
var = sigma**2
# Managing the language
if lang == 'en':
    labelWigner = r'Wigner semicircle distribution, excess kurtosis $=-1$'
    labelNormal = r'Normal distribution, excess kurtosis $=0$'
    labelHypSec = r'Hyperbolic secant distribution, excess kurtosis $=2$'
    title = ('Right tails of various symmetric distributions\n'
              + r'They all have the same standard deviation ($\sigma=1$)')
elif lang == 'fr':
    labelWigner = r'Loi du demi-cercle de Wigner, excès de kurtosis $=-1$'
    labelNormal = r'Loi normale, excès de kurtosis $=0$'
    labelHypSec = r'Loi sécante hyperbolique, excès de kurtosis $=2$'
    title = ''
# Create a new plot
fig, ax = plt.subplots(figsize = set_size())
# Plot all these distributions
ax.plot(x, wigner(x, sigma), label = labelWigner)
ax.plot(x, gaussian(x, mu, sigma), label = labelNormal)
ax.plot(x, sech(x), label = labelHypSec)
# Graph glitter
ax.set(xlabel=r'$x$', ylabel=r'$f(x)$')
if title:
    ax.set(title = title)
ax.grid()
ax.set_ylim(bottom = 0)
ax.set_xlim(left = min(x), right = max(x))
plt.legend(framealpha = 1.0)
fig.savefig("Kurtosis.pdf", bbox_inches='tight')
# Clean-up for next plot
plt.close('all')

# Focus on ]0;2]
x = np.linspace(0, 2, 1000)[1:]
# Managing the language
if lang == 'en':
    label = 'Skewness'
    title=('Log-normal distributions for various skewnesses (solid lines) '
           + 'and\ntheir normal counterpart (dashed lines, same standard '
           + 'deviation)\n'+ r'They all have the same median$= 1$')
elif lang == 'fr':
    label = 'Asymétrie'
    title = ''
# Create a new plot
fig, ax = plt.subplots(figsize = set_size())
# Find sigmas corresponding to our skewnesses (a cubic problem)
skewnesses = np.array([0.5, 1.0, 1.5])
a = 1
b = 3
c = 0
d = -4 - skewnesses**2
D = 18*a*b*c*d - 4*b**3*d + b**2*c**2 - 4*a*c**3 - 27*a**2*d**2
if not np.all(D < 0):
    raise Exception('A log-normal distribution can only have 1 skewness')
D0 = b**2 - 3*a*c
D1 = 2*b**3 - 9*a*b*c + 27*a**2*d
C = np.cbrt((D1 + np.sqrt(D1**2 - 4*D0**3))/2)
sol = -(b + C + D0/C)/(3*a)
sigmas = np.sqrt(np.log(sol))
# Plot it
for sigma, skewness in zip(sigmas, skewnesses):
    color = next(ax._get_lines.prop_cycler)['color']
    ax.plot(x, gaussian(np.log(x), mu, sigma)/x, color = color,
            label = label + ' = ' + str(skewness))
    ax.plot(x, gaussian(x, 1, sigma), color = color, linestyle = '--',
            linewidth = 1.0)
# Graph glitter
ax.set(xlabel=r'$x$', ylabel=r'$f(x)$')
if title:
    ax.set(title = title)
ax.grid()
ax.set_ylim(bottom = 0)
ax.set_xlim(left = min(x), right = max(x))
plt.legend(framealpha = 1.0)
fig.savefig("Lognormal.pdf", bbox_inches='tight')
# Clean-up for next plot
plt.close('all')
