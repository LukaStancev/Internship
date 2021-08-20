#
#  Basic functions, useful for a wide variety of plots
#  Usage : from BasicFunctions import *
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import moment
import serpentTools
from serpentTools.settings import rc

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

# Read CSV file with comments
# Source : https://stackoverflow.com/questions/14158868/python-skip-comment-lines-marked-with-in-csv-dictreader
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
#  Functions for proper formatting of Matpotlib figures in LaTeX
#  See: https://jwalton.info/Embed-Publication-Matplotlib-Latex/
#---
def set_size(aspect = 'default'):
    # Output of \showthe\textwidth in LaTeX
    width = 472.03111 # in points
    width = width / 72.27 # into inches
    if aspect == 'halfsquare':
        width = width/2
    if aspect == 'default':
        # Retrieve default figsize
        default = plt.rcParams["figure.figsize"]
        ratio = default[1] / default[0]
    elif aspect == 'square' or aspect == 'halfsquare':
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
