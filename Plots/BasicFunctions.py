#
#  Basic functions, useful for a wide variety of plots
#  Usage : from BasicFunctions import *
#  Author ; V. Salino (IRSN), 02/2021
#

# Imports
import numpy as np

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
