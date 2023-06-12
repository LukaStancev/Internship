#
#  Basic functions, useful for a wide variety of plots
#  Usage : from BasicFunctions import *
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import lcm
import os
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy.stats import moment
#import serpentTools
#from serpentTools.settings import rc

# Layout of the assemblies in the core:
#   = 0, void (i.e. no assembly at this position)
#   = 1, presence of an assembly
def GetCoreLayout(nbassembly):
    if nbassembly == 157: # assemblies
        CoreLayout = np.array([[1, 1, 0, 0, 0, 0, 0, 0],
                               [1, 1, 1, 1, 0, 0, 0, 0],
                               [1, 1, 1, 1, 1, 0, 0, 0],
                               [1, 1, 1, 1, 1, 1, 0, 0],
                               [1, 1, 1, 1, 1, 1, 1, 0],
                               [1, 1, 1, 1, 1, 1, 1, 0],
                               [1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1]], dtype=int)
    else:
        raise Exception('Unsupported core layout')
    FullCoreLayout = UnfoldQuarter(CoreLayout)
    return CoreLayout, FullCoreLayout

# Layout for Serpent output tallies
# outside : number of assembly-sized tallies outside the core, on each side
def GetSerpentLayout(FullCoreLayout, outside):
    hstack = np.zeros((len(FullCoreLayout), outside))
    SerpentLayout = np.hstack([hstack, FullCoreLayout, hstack])
    vstack = np.zeros((outside, len(FullCoreLayout) + 2*outside))
    SerpentLayout = np.vstack([vstack, SerpentLayout, vstack])
    return SerpentLayout

# Retrieve a Serpent power distribution
def ReadPdistrSerpent(file, FullCoreLayout, msg = '', reaction = '-80'):
    # Retrieve Serpent detector output
    det = serpentTools.read(file)
    # Default (-80) is total energy deposition, see
    # http://serpent.vtt.fi/mediawiki/index.php/ENDF_reaction_MT%27s_and_macroscopic_reaction_numbers
    power = det.detectors[reaction].tallies
    if reaction == '-80':
        power = power[0]
    # Number of assembly-sized tallies outside the core, on each side
    outside = 2
    # Produce the core layout corresponding to tallies
    hstack = np.zeros((len(FullCoreLayout), outside))
    TalliesLayout = np.hstack([hstack, FullCoreLayout, hstack])
    vstack = np.zeros((outside, len(FullCoreLayout) + 2*outside))
    TalliesLayout = np.vstack([vstack, TalliesLayout, vstack])
    # Compute fraction of power delivered in the reflector
    corepower = np.where(TalliesLayout == 0, 0, power[1])
    radialreflpower = np.where(TalliesLayout == 1, 0, power[1])
    botreflpower = np.sum(power[0])
    topreflpower = np.sum(power[2])
    reflpower = np.sum(radialreflpower) + botreflpower + topreflpower
    print('Fraction of the power Serpent delivers in the reflector = '
          + str(np.sum(reflpower)/np.sum(power)) + ' (' + msg + ')')
    # Ignore reflector power and renormalize core power
    corepower = corepower/np.sum(corepower)*np.count_nonzero(corepower)
    # The input is symmetrical by eighth, so we can symmetrize (a fold followed
    # by an unfold) in order to increase the statistical strength
    symcorepower = UnfoldEighth(FoldEighth(corepower))
    # Estimate the uncertainty
    sym = symcorepower
    unsym = corepower
    reldif = unsym[np.nonzero(unsym)] / sym[np.nonzero(sym)] - 1
    print("min, max (%)=" + str(np.min(reldif)*100) + ", "
          + str(np.max(reldif)*100))
    # Return powers, removing all of them equal to zero
    return symcorepower[np.nonzero(symcorepower)]

# Retrieve a per-batch Serpent power distribution
def ReadPdistrHistSerpent(file, FullCoreLayout, reaction = '-80'):
    SerpentLayout = GetSerpentLayout(FullCoreLayout, 2)
    # Read that Serpent history file
    hist = serpentTools.read(file)
    # For a detector (here -80, by default), the results are grouped in three
    # columns that provide (1) the cycle-wise value,
    #                      (2) the cumulative mean and
    #                      (3) the corresponding relative statistical error.
    # Cf. http://serpent.vtt.fi/mediawiki/index.php/Description_of_output_files#History_output
    # Only the first is kept here.
    cycle = hist['det' + reaction][:, 0::3]
    if reaction == '-80':
        # Among that, only the total deposited energy is kept
        cycle = cycle[:, 0::3]
    # Among that, only the central layer (the active core) is kept ; top and
    # bottom layers are splitted out
    cycle = np.split(cycle, 3, axis = 1)[1]
    # Raise exception if inactive generations are being used. One (1) inactive
    # generation is the minimum permitted by Serpent and therefore the only
    # one we accept.
    if np.isclose(0.0, np.sum(cycle[1]), atol = 1e-20):
        raise Exception('The convergence plot cannot be achieved since '
                        + 'inactive generations have been detected: '
                        + '*set pop [x] [y] 1* should be used in Serpent.')
    # Reshaping into nbatch*size*size
    nbatch = np.shape(cycle)[0]
    nxysize = int(math.sqrt(np.shape(cycle)[1]))
    cycle.shape = (nbatch, nxysize, nxysize)
    # Prepare output
    out = np.zeros((nbatch, np.count_nonzero(FullCoreLayout)))
    for ibatch in range(1, nbatch):
        # Filter out the reflector to get only the response of the core
        thiscycle = np.where(SerpentLayout == 0, 0, cycle[ibatch])
        # Normalize the core power
        nassembly = np.count_nonzero(FullCoreLayout)
        thiscycle = nassembly * thiscycle / np.sum(thiscycle)
        # Remove columns and rows full of zeroes
        thiscycle = Deploy2D(thiscycle[np.nonzero(thiscycle)], FullCoreLayout)
        # The input is symmetrical by eighth, so we can symmetrize (a fold
        # followed by an unfold) in order to increase the statistical strength
        thiscycle = UnfoldEighth(FoldEighth(thiscycle))
        # Keep all the (non-zero) powers of the full core
        out[ibatch] = thiscycle[np.nonzero(thiscycle)]
        # ... except the first (empty) one (it is the only inactive generation)
    return out[1:]

# Unfold a quarter core on a full core layout
def UnfoldQuarter(quarter):
    # Prepare upper core layout, minus central vertical line
    UpperLeft = np.delete( np.rot90(quarter), -1, axis=0)
    UpperRight = np.delete( np.delete(quarter, -1, axis=0), 0, axis=1 )
    Upper = np.hstack( (UpperLeft, UpperRight) )
    # Prepare lower core layout, letting vertical line
    LowerLeft = np.rot90(quarter, 2)
    LowerRight = np.delete( np.rot90(quarter, 3), 0, axis=1 )
    Lower = np.hstack( (LowerLeft, LowerRight) )
    return(np.vstack((Upper, Lower)))

# Fold in quarter maps (averaging per rotation) for plotting purposes
def fold(full):
    # Retrieve center coordinates (assuming a square)
    center = len(full[:, 0])//2
    # We retrieve each quarter, that are not squares (+1 in one dimension)
    UpperRight = full[:center, center:]
    LowerRight = np.rot90(full[center:, center+1:])
    LowerLeft = np.rot90( full[center+1:, :center+1], 2)
    UpperLeft = np.rot90( full[:center+1, :center], 3)
    # Cumulate every quarters and mean them
    Quarter = (UpperRight + LowerRight + LowerLeft + UpperLeft)/4
    # Add the horizontal line. The information contained here is redundant,
    # except for the central value (which is full[center, center])
    LowerLine = np.hstack(( full[center, center], np.flip(Quarter[:, 0]) ))
    return np.vstack(( Quarter, LowerLine ))

# Unfold a eighth map into a full map
def UnfoldEighth(eighth):
    # Complete the eighth so that it forms a square, the future quarter
    nrows = eighth.shape[0]
    ncols = eighth.shape[1]
    if nrows != ncols:
        eighth = np.hstack((eighth, np.zeros((nrows, nrows - ncols))))
    # Produce a quarter map (southeast)
    quarter = np.zeros((nrows, nrows))
    for x in range(0, nrows):
        for y in range(0, nrows):
            if x < y:
                quarter[x, y] = eighth[y, x]
            else:
                quarter[x, y] = eighth[x, y]
    # Rotate southeast quarter into northeast quarter
    quarter = np.rot90(quarter)
    # Unfold to full map
    full = UnfoldQuarter(quarter)
    return(full)

# Fold a full map into an eighth
def FoldEighth(full):
    quarter = np.rot90(fold(full), 3)
    size = quarter.shape[0]
    for x in range(1, size):
        for y in range(1, x):
            quarter[x, y] = (quarter[x, y] + quarter[y, x])/2
            quarter[y, x] = 0
    return(quarter)

# Deploy a 1D-indexed map into a 2D-indexed map, on a given layout
def Deploy2D(OneD, FullCoreLayout):
    TwoD = np.zeros(( len(FullCoreLayout[:, 0]),   # Y-axis
                      len(FullCoreLayout[0, :]) )) # X-axis
    i = 0
    for x in range(0, len(FullCoreLayout[0, :])):
        for y in range(0, len(FullCoreLayout[:, 0])):
            if FullCoreLayout[x, y] != 0:
                TwoD[x, y] = OneD[i]
                i = i + 1
    return TwoD

# Sum over submatrices of size n x n.
def matrixaverage(matrix, n):
    # https://moonbooks.org/Articles/How-to-downsampling-a-matrix-by-averaging-elements-nn-with-numpy-in-python-/Edit/
    b = matrix.shape[0]//n
    return matrix.reshape(-1, n, b, n).sum((-1, -3))/n

#---
# Read CSV file with comments
# Source : https://stackoverflow.com/questions/14158868/python-skip-comment-lines-marked-with-in-csv-dictreader
#---
def decomment(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw

#---
#  Functions useful for calculating standard errors (SE) for standard
#  deviation (SD), skewness (SK) and excess kurtosis (KU)
#---
# https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
def SE_Variance(sample, std):
    n = len(sample)
    fourth = moment(sample, moment = 4)
    return np.sqrt(1/n*(fourth - (n - 3)/(n - 1)*std**4))
def SE_SD(sample):
    std = np.std(sample, ddof = 1)
    return SE_Variance(sample, std)/(2*std)
# https://www.researchgate.net/publication/285590690_Standard_errors_A_review_and_evaluation_of_standard_error_estimators_using_Monte_Carlo_simulations
def SE_SD_Normal(n):
    return 1/np.sqrt(2*(n-1))
def SE_SK(n):
    return np.sqrt( (6*n*(n-1)) / ((n - 2)*(n + 1)*(n + 3)) )
def SE_KU(n):
    return 2*SE_SK(n)*np.sqrt( (n**2 - 1) / ((n - 3)*(n + 5)) )

#---
#  Equivalent to NumPy percentile (or quantile), but supports weights
#  Ref. : https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
#---
def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)
#---
#  Compute weighted skewness and excess kurtosis
#  Ref. : https://stackoverflow.com/questions/61521371/calculate-weighted-statistical-moments-in-python
#---
def weighted_mean(var, wts):
    """Calculates the weighted mean"""
    return np.average(var, weights=wts)

def weighted_variance(var, wts):
    """Calculates the weighted variance"""
    return np.average((var - weighted_mean(var, wts))**2, weights=wts)

def weighted_std(var, wts):
    """Calculates the weighted standard deviation"""
    return np.sqrt(weighted_variance(var, wts))

def weighted_relstd(var, wts):
    """Calculates the weighted relative standard deviation"""
    return weighted_std(var, wts)/weighted_mean(var, wts)*100

def weighted_skew(var, wts):
    """Calculates the weighted skewness"""
    return (np.average((var - weighted_mean(var, wts))**3, weights=wts) /
            weighted_variance(var, wts)**(1.5))

def weighted_kurtosis(var, wts):
    """Calculates the weighted skewness"""
    return (np.average((var - weighted_mean(var, wts))**4, weights=wts) /
            weighted_variance(var, wts)**(2) - 3)
#---
#  Compute weighted correlation matrix
#  Ref. : https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas
#         https://numpy.org/doc/stable/reference/generated/numpy.cov.html
#---
def weighted_corrcoef(weights, x, y = None):
    """Weighted Correlation"""
    if x.ndim == 1:
        # Only two samples (x and y)
        corr = (np.cov(x, y, aweights = weights) /
                np.sqrt(np.cov(x, x, aweights = weights) *
                np.cov(y, y, aweights = weights)))
    else:
        # N samples, all stored in x
        cov = np.cov(x, aweights = weights)
        stddev = np.sqrt(np.diag(cov))
        stddevproduct = np.outer(stddev, stddev)
        corr = cov/stddevproduct
    return corr

#---
#  Round a correlation matrix
#---
def roundcorr(correlation, Energies, delta):
    ngrp = len(Energies) - 1
    correlation = np.around(correlation/delta)*delta
    # We may have some data without variation, such as
    # * the light water TSLs,
    # * the inelastic O16 and B10 cross sections, or
    # * the threshold reactions.
    # The correlation cannot be computed in such cases because their standard
    # deviation (in the correlation's denominator) is zero. We replace these
    # NaN with zeros.
    correlation = np.nan_to_num(correlation, nan = 0)
    # Sequentialize the correlations (flatten matrix into vector)
    seqcor = np.ravel(correlation, order = 'C')
    xlim = np.unique((np.where(seqcor[:-1] != seqcor[1:])[0] + 1) % ngrp)
    # Now that this has been done on the columns, do the same thing on the rows
    seqcor = np.ravel(correlation, order = 'F')
    ylim = np.unique((np.where(seqcor[:-1] != seqcor[1:])[0] + 1) % ngrp)
    # If the matrix is totally homogeneneous (a very exceptionnal yet possible
    # case), add the zero index
    if 0 not in xlim:
        xlim = np.concatenate(([0], xlim))
    if 0 not in ylim:
        ylim = np.concatenate(([0], ylim))
    # Now that we've found the homogeneous rectangles, we keep only one value
    # (the first here, since they are all equal)
    XX, YY = np.meshgrid(xlim, ylim)
    correlation = correlation[YY, XX]
    # For the energies, add the final point (energy mesh has G+1 points)
    xEnergies = np.concatenate((Energies[xlim], [Energies[-1]]))
    yEnergies = np.concatenate((Energies[ylim], [Energies[-1]]))
    return correlation, xEnergies, yEnergies

#---
#  Functions for proper formatting of Matpotlib figures in LaTeX
#  See: https://jwalton.info/Embed-Publication-Matplotlib-Latex/
#---
def set_size(aspect = 'default', bonus = False):
    # Output of \showthe\textwidth in LaTeX
    width = 472.03111 # in points
    width = width / 72.27 # into inches
    if aspect == 'halfsquare':
        width = width/2
    if aspect == 'default':
        # Retrieve default figsize in order to keep default proportionnality
        default = plt.rcParams["figure.figsize"]
        ratio = default[1] / default[0]
    elif aspect == 'square' or aspect == 'halfsquare':
        ratio = 1
        # An incompatibility between savefig's bbox_inches='tight' and
        # subplots' figsize forces us to use dirty fixed ratios.
        # See :
        # https://github.com/matplotlib/matplotlib/issues/11681
        # https://stackoverflow.com/questions/16118291/matplotlib-make-final-figure-dimensions-match-figsize-with-savefig-and-bbox-e
        # https://stackoverflow.com/questions/16032389/pad-inches-0-and-bbox-inches-tight-makes-the-plot-smaller-than-declared-figsiz
        if aspect == 'square':
            if bonus:
                #width = width*16.59/14.06*58.7/59.1*39.35/39.95
                width = width*1.1543556
            else:
                #width = width*16.59/13.48
                width = width*1.2307121
        elif aspect == 'halfsquare':
            #width = width*29.3/26.45*1.015
            width = width*1.1243667
    elif aspect == 'fullA4':
        #width = width*23.3/18.1*24.15/28.025
        width = width*1.1092996
        default = plt.rcParams["figure.figsize"]
        #ratio = default[1] / default[0] / 1.1092996 * (20-2)/8.4 * 18/22.2
        ratio = default[1] / default[0] * 1.566260131
    elif aspect == 'fullA4-Draglib':
        width = width*1.0785539
        default = plt.rcParams["figure.figsize"]
        ratio = default[1] / default[0] * 1.649068277
    elif aspect == 'slide-Draglib':
        width = width*1.0785539
        default = plt.rcParams["figure.figsize"]
        ratio = default[1] / default[0] * 0.959847297
        #ratio = default[1] / default[0] * 1.649068277 *(14.5/17.5)/(24.2/17)
    elif aspect == 'speedup':
        #width = width*24.45/25.6
        width = width*1.17542620488
        ratio = 1
    elif aspect == 'correlation':
        #width = width*23.3/22.3*22.3/25.8
        width = width*0.90310077
        ratio = 1
    else:
        raise Exception('"' + aspect + '" aspect has not been defined.')
    # Keep the default proportionality
    height = width * ratio
    return (width, height)
def tex_fonts():
    tex_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 12pt font in plots, to match 12pt font in document
        "axes.labelsize": 12,
        "font.size": 12,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    }
    return tex_fonts

#---
#  Plotting functions
#---
def FullCoreGlitter(ax, FullCoreLayout):
    # Remove spines
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # Add squares around each assembly
    DrawSquares(FullCoreLayout, ax)
    # Add Battleship-style coordinates as ticks
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    xlabels = ['R', 'P', 'N', 'M', 'L', 'K', 'J', 'H', 'G', 'F', 'E', 'D', 'C',
               'B', 'A']
    ylabels = list(np.arange(1, 16))
    plt.setp(ax, xticks = np.arange(15), xticklabels = xlabels,
             yticks = np.arange(15), yticklabels = ylabels)
    ax.tick_params(axis = 'both', which = 'both', length = 5, color = 'w')
# Add lines around each square
def DrawSquares(FullCoreLayout, ax):
    i = 0
    for line in FullCoreLayout:
        xmin, xmax = LayoutProportions(FullCoreLayout, i)
        if i > 0:
            xmin_above, xmax_above = LayoutProportions(FullCoreLayout, i - 1)
            xmin = min(xmin, xmin_above)
            xmax = max(xmax, xmax_above)
        ax.axhline(y = -0.5 + i, xmin = xmin, xmax = xmax, color = 'black',
                   linewidth = 1.0)
        ax.axvline(x = -0.5 + i, ymin = xmin, ymax = xmax, color = 'black',
                   linewidth = 1.0)
        i = i + 1
    # Add the last bottom line
    xmin, xmax = LayoutProportions(FullCoreLayout, i - 1)
    ax.axhline(y = -0.5 + i, xmin = xmin, xmax = xmax, color = 'black',
               linewidth = 1.0)
    ax.axvline(x = -0.5 + i, ymin = xmin, ymax = xmax, color = 'black',
               linewidth = 1.0)
def LayoutProportions(FullCoreLayout, i):
    line = FullCoreLayout[i]
    # Count the number of zeroes on each side
    zeroes = [sum(g) for j, g in itertools.groupby(line+1) if j==1]
    if zeroes == []:
        zeroes = [0, 0]
    # Count the ones in the center
    ones = sum(line)
    total = sum(zeroes, ones)
    if total != len(FullCoreLayout):
        raise Exception('Review the line number ' + str(i) +' of core layout. '
                        + 'It is not as regular as expected : ' + str(line))
    xmin = zeroes[0]/total
    xmax = (zeroes[0] + ones)/total
    return (xmin, xmax)
# Add visual aids, such as median and diagonal crosses
def VisualAids(ax):
    # Common plotting options
    lwidth = 1.0
    col = "black"
    alpha = 0.3
    # Medians
    ax.plot([0.5, 0.5], [0, 1], transform=ax.transAxes,
            linewidth = lwidth, color = col, alpha = alpha)
    ax.plot([0, 1], [0.5, 0.5], transform=ax.transAxes,
            linewidth = lwidth, color = col, alpha = alpha)
    # Diagonals
    ax.plot([0, 1], [0, 1], transform=ax.transAxes,
            linewidth = lwidth, color = col, alpha = alpha)
    ax.plot([1, 0], [0, 1], transform=ax.transAxes,
            linewidth = lwidth, color = col, alpha = alpha)
# https://stackoverflow.com/questions/45441909/how-to-add-a-fixed-width-border-to-subplot
def add_subplot_border(ax, width=1, color=None, zorder=1):

    fig = ax.get_figure()

    # Convert bottom-left and top-right to display coordinates
    x0, y0 = ax.transAxes.transform((0, 0))
    x1, y1 = ax.transAxes.transform((1, 1))

    # Convert back to Axes coordinates
    x0, y0 = ax.transAxes.inverted().transform((x0, y0))
    x1, y1 = ax.transAxes.inverted().transform((x1, y1))

    rect = plt.Rectangle(
        (x0, y0), x1-x0, y1-y0,
        color=color,
        zorder=zorder,
        transform=ax.transAxes,
        lw=2*width+1,
        fill=None,
    )
    fig.patches.append(rect)
#---
#  lcm loading function, dealing with PyGan's peculiarities
#---
def loadlcm(filepath):
    if '/' in filepath:
        # PyGan can only read files in the local directory. When the file to be
        # read is located elsewhere, it is necessary to make a link.
        filename = filepath.split('/')[-1]
        os.system('ln -s ' + filepath + ' _' + filename)
        # PyGan's lcm method requires that filenames begin with an underscore,
        # but at the same time requires that it is passed without this
        # underscore
        loadedfile = lcm.new('LCM_INP', filename)
        # Once loaded, we can delete the created link
        os.system('rm _' + filename)
    else:
        loadedfile = lcm.new('LCM_INP', filepath)
    return loadedfile
