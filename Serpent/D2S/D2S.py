#  Convert library-type and geometry-type DRAGON objects to Serpent input file.
#  It is limited to PWR geometries symmetrical by eighth.
#  Usage   : from D2S import D2S
#            D2S(geofilename)
#  Authors : V. Salino, B. Dechenaux Sensei
#  Date    : 11/2020

#---
# Imports
#---
import math
import numpy as np
import os
import decimal
# faulthandler for segmentation faults. Left here due to a heisenbug in lcm.
import faulthandler; faulthandler.enable()
import lcm

#---
#  Classes
#---
class compo:
    """
    The compo class stores the isotopic composition of a material.
    """
    def __init__(self, filename, mix, temp, material):
        # Original filename from which mix was extracted. Several files often
        # contain identical compositions.
        self.m_filename = filename
        # Dictionary containing the isotopes (key) and their atomic densities
        # in at/(b.cm) (values)
        self.m_compo = material
        # Integer designating the mix number
        self.m_mix = mix
        # Temperature in Kelvin - applicable to all isotopes
        if len(set(temp)) != 1:
            raise Exception('MIX ' + str(mix) + ' in ' + filename +
                    ' has several temperatures. Serpent has a limitation of ' +
                    'only one temperature per MIX.')
        self.m_temp = temp[0]
        # Forget about negligible concentrations at initialization
        self.trimCompo(1.0E-9)
    def getCompo(self):
        return self.m_compo
    def getMix(self):
        return self.m_mix
    def getTemp(self):
        return self.m_temp
    def getFilename(self):
        return self.m_filename
    def trimCompo(self, threshold):
        """
        Modifies the dictionary of compo values to cut below a threshold value
        """
        newCompo = {}
        for key,value in self.m_compo.items():
            if value > threshold:
                newCompo[ key ] = value
        self.m_compo = newCompo
    def writeSerpent(self, out):
        """
        Writes the isotopic composition in a Serpent input file
        """
        # Declare S(a,b) for H1_H2O, if it present. Here, Serpent requires the
        # ZAID :
        # http://serpent.vtt.fi/mediawiki/index.php/Input_syntax_manual#mat_moder
        tslFiles = None
        moder = ''
        for iso in self.m_compo:
            if iso == 'H1_H2O':
                moder = 'moder lwtr' + str(self.m_mix) + ' 1001 '
        # Write the header line for the entire material...
        out.write('mat mix' + str(self.m_mix) + ' sum ' + moder
                  + 'tmp ' + str(self.m_temp) + ' % Kelvin\n')
        # ...and then, write one line for each isotope. But let's prepare it,
        # first !
        for iso in self.m_compo:
            if iso == 'H1_H2O':
                iso_ace = 'H1lwtr'
                iso_ace_tsl = 'lwtr'
            else:
                iso_ace = iso
            # Subset xsdata file for this isotope we're on, for all available
            # temperatures
            xsdata_subset = []
            # Go through xsdata file, find ace files of this isotope we're on
            with open('../../../Njoy/Universal.xsdata') as xsdatafile:
                xsdatalines = xsdatafile.readlines()
                for xsdataline in xsdatalines:
                    # First field in xsdata file contains isotope's name
                    isoxsdata = xsdataline.split()[0]
                    if isoxsdata.startswith(iso_ace + '.'):
                        xsdata_subset.append(xsdataline.rstrip().split())
            if len(xsdata_subset) == 0:
                raise Exception('Could not find ace file for isotope: ' + iso)
            # Sort by ascending temperatures, contained in field [6]
            xsdata_subset.sort(key=lambda x: int(x[6]))
            # Check that a temperature below the one requested is available
            if int(xsdata_subset[0][6]) > self.m_temp:
                raise Exception('The minimum temperature available in the ace '
                + 'files of the isotope ' + iso_ace + ' is '
                + xsdata_subset[0][6] + 'K' + '. This is too high for the '
                + 'requested temperature: ' + str(self.m_temp) + 'K')
            # For TSL interpolation, we also need an available temperature
            # above the one requested
            if (int(xsdata_subset[-1][6]) < self.m_temp) and (iso == 'H1_H2O'):
                raise Exception('The maximal temperature available in the TSL '
                + 'ace files of the isotope ' + iso_ace + ' is '
                + xsdata_subset[-1][6] + 'K' + '. This is too low for the '
                + 'requested temperature: ' + str(self.m_temp) + 'K')
            # Select ace file with the temperature immediately below
            aceFile = None
            i = 0
            while not aceFile:
                if int(xsdata_subset[i][6]) > self.m_temp:
                    aceFile = xsdata_subset[i-1][1]
                    # Corresponding TSL files should also be kept (lower *and*
                    # upper bounds are required)
                    if iso == 'H1_H2O':
                        if tslFiles:
                            raise Exception('TSL file has been already '
                            + 'attributed. D2S is limited to one single TSL '
                            + 'per material.')
                        # * Remove isotope name used in continuous ace and
                        #   use the TSL specific name instead
                        # * Remove last character ('c' for 'continuous ace')
                        #   and replace it with 't' (for 'tsl ace')
                        tslFiles = (iso_ace_tsl + xsdata_subset[i-1][1][-4:-1]
                                    + 't' + ' '
                                    + iso_ace_tsl + xsdata_subset[i][1][-4:-1]
                                    + 't')
                i = i + 1
            # Write the line for the isotope we're on
            out.write("%s %.8E\n" %(aceFile,self.m_compo[iso]))
        # Write the 'thermr' card, if needed
        if tslFiles:
            out.write('therm lwtr' + str(self.m_mix) + ' ' + str(self.m_temp)
                      + ' ' + tslFiles + '\n')
        out.write('\n')
    def __eq__(self, compoToCompare):
        """
        Overloading == operator : we can compare 2 objects directly
        """
        if self.m_compo == compoToCompare.getCompo():
            return True
        else:
            return False
    def __repr__(self):
        """
        Display a composition mix and original file
        """
        return "compo of mix %d from file %s"%(self.m_mix,self.m_filename)
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
#---
#  Math functions
#---
# Perform a truncation, keeping only n significant digits
def truncate(x, n):
    def fexp(number):
        (sign, digits, exponent) = decimal.Decimal(number).as_tuple()
        return len(digits) + exponent - 1
    return math.floor(x/(10**(fexp(x)-n+1)))*(10**(fexp(x)-n+1))
# These exact trigonometric functions avoid the use of Pi and therefore
# rounding off errors such as sin(pi) = 1.0e-16
def exact_sin(angle):
    if angle % 180 == 0:
        result = 0
    elif angle % 360 == 90:
        result = 1
    elif angle % 360 == 270:
        result = -1
    else:
        raise Exception('exact_sin function was not designed for angle=',angle)
    return result
def exact_cos(angle):
    if angle % 360 == 0:
        result = 1
    elif angle % 180 == 90:
        result = 0
    elif angle % 360 == 180:
        result = -1
    else:
        raise Exception('exact_cos function was not designed for angle=',angle)
    return result
#---
#  Composition function
#---
def getMat(idx, mix, dens, name):
    if len(name) != len(dens):
        raise Exception('ISOTOPERNAME and ISOTOPESDENS have different lengths '
        + 'but should not. It is most probably due to ISOTOPERNAME split.')
    mat = {}
    ind = np.where(mix == idx)
    # Loop over mixes, stored in this list: ind[0]
    for i in ind[0]:
        mat[ name[i] ] = dens[i]
    return mat
#---
#  Geometric functions
#---
# Compute number of pins and performs checks regarding that matter
def getNpin(npinx, npiny, cellIDs):
    if npinx != npiny:
        raise Exception('Non-square lattices are not supported.')
    # We're now sure it is a square so we can keep only one dimension
    npin = npinx
    # Checking number of pins with Gauss formula
    if len(cellIDs) != npin*(npin+1)/2:
        raise Exception('Non-eighth geometries are not supported.')
    # Compute number of pins on the full pin map (on each axis)
    npin = (npin-1)*2+1
    return npin
# Deploy cell numbers on the full pin map
def deploy(npin, cellIDs):
    center = int((npin-1)/2)
    center = np.array([center, center])
    # Initialize pin map
    pinmap = np.zeros((npin,npin))
    i = 0
    # cellIDs are expressed on east-south-east eighth, that is
    # 0 0 0 0 0
    # 0 0 0 0 0
    # 0 0 x x x
    # 0 0 0 x x
    # 0 0 0 0 x
    for x in range(int((npin-1)/2), npin):
        for y in range(x, npin):
            # cellIDs are negative because they refer to embedded geometries
            # We keep the positive integers instead
            pinmap[x][y] = -cellIDs[i]
            i = i + 1
            xp = x
            yp = y
            for angle in range(0, 360-45, 45):
                # https://en.wikipedia.org/wiki/Rotations_and_reflections_in_two_dimensions
                # For each angle, we calculate the matrix that will operate the
                # reflection relative to a line forming an angle with the
                # horizontal, in order to reconstruct every eigth in the
                # counter-clockwise order
                Refl = np.array([
                    [exact_cos(2*angle), exact_sin(2*angle)],
                    [exact_sin(2*angle),-exact_cos(2*angle)]])
                # Reflection line goes through the assembly center
                [xp,yp] = np.dot(np.array([xp,yp]) - center, Refl) + center
                pinmap[xp][yp] = pinmap[x][y]
    return pinmap
# Compute pin and assembly pitches
def getPitches(npin, pinmap, cellNames, meshxList, meshyList):
    # Pick first pincell that is not a border pincell, i.e. pin[1][1] and
    # retrieve its pitch as a reference against which the other pitches will be
    # checked
    cellID = int(pinmap[1][1])
    cellName = cellNames[cellID-1]
    pinpitch = meshxList[cellName]
    assemblypitch = 0
    # Loop over every pincells
    for [x,y], cellID in np.ndenumerate(pinmap):
        cellID = int(cellID)
        cellName = cellNames[cellID-1]
        meshx = meshxList[cellName]
        meshy = meshyList[cellName]
        # Compute assembly pitch on first column
        if y == 0:
            assemblypitch = assemblypitch + meshy
        # Perform consistency checks on X and Y pincell sizes
        if (0 < x < npin-1) and (0 < y < npin-1):
            # Central pincells
            if not math.isclose(meshx, pinpitch) or \
               not math.isclose(meshy, pinpitch):
                raise Exception('Irregular lattices are not supported.')
        else:
            if (0 < x < npin-1) or (0 < y < npin-1):
                # Border pincells that we allow to have a slighly larger size,
                # in order to allow a water gap
                if not (pinpitch*(1-1e-6) < meshx < pinpitch*1.1):
                    raise Exception('Irregular lattices are not supported.')
                if not math.isclose(meshy, pinpitch):
                    raise Exception('Irregular lattices are not supported.')
            else:
                # Corner pincells that we allow to have a slighly larger size
                # in both directions, in order to allow a water gap
                if not (pinpitch*(1-1e-6) < meshx < pinpitch*1.1):
                    raise Exception('Irregular lattices are not supported.')
                if not (pinpitch*(1-1e-6) < meshy < pinpitch*1.1):
                    raise Exception('Irregular lattices are not supported.')
    return pinpitch,assemblypitch
# Main function
def D2S(filepath):
    global compo
    filename = filepath.split('/')[-1]
    #---
    #  Prepare file writing for Serpent
    #---
    out = open(filename + '.sss2', 'w')
    #---
    #  Write general information
    #---
    out.write('set title ' + '"Tihange ' + filename + '"\n')
    out.write('set acelib "../../Njoy/Universal.xsdata"\n')
    out.write('set bc 3\n')
    out.write('set pop 6000 500 20\n')
    out.write('plot 3 2500 2500\n')
    out.write('mesh 3 2500 2500\n\n')
    #---
    #  Retrieve assembly geometry
    #---
    fuel_geo = loadlcm(filepath + '.geo')
    # Retrieve name given to each type of cell
    cellNames = fuel_geo['CELL'].split()
    # Retrieve identifier of each pincell (i.e. the assembly layout)
    cellIDs = fuel_geo['MIX']
    # Retrieve number of pins in both X-axis and Y-axis (on eighth geometry)
    npinx = fuel_geo['STATE-VECTOR'][2]
    npiny = fuel_geo['STATE-VECTOR'][3]
    # Retrieve pincell sizes
    meshxList = {}
    meshyList = {}
    for cellName in cellNames:
        meshxList[cellName] = fuel_geo[cellName]['MESHX'][1]
        meshyList[cellName] = fuel_geo[cellName]['MESHY'][1]
    # Perform consistency checks regarding number of pins
    npin = getNpin(npinx,npiny,cellIDs)
    # Deploy the full pin map from eighth vectorized cellIDs
    pinmap = deploy(npin,cellIDs)
    # Compute pin and assembly pitches in a thorough consistency check
    [pinpitch, assemblypitch] = getPitches(npin,pinmap,cellNames,meshxList,
                                           meshyList)
    # Assess the compositions used in this geometry
    mix_in_geom = set()
    for cellID in range(0,len(cellNames)):
        cellName = cellNames[cellID]
        materials = fuel_geo[cellName]['MIX']
        mix_in_geom.update(materials.tolist())
    mix_in_geom = list(mix_in_geom)
    mix_in_geom.sort()
    #---
    #  Write assembly geometry
    #---
    # Writing the few different type of pincells
    for cellID in range(0,len(cellNames)):
        cellName = cellNames[cellID]
        radii = fuel_geo[cellName]['RADIUS'][1:]
        materials = fuel_geo[cellName]['MIX']
        # Writes cellID together with its name as a commentary
        out.write('pin ' + str(cellID+1) + ' %' + cellName + '\n')
        for radius,material in zip(radii, materials):
            # I would prefer to align on the largest material name, in
            # particular for debugging purposes, but... How to?
            out.write('mix' + str(material) + ' ' + str(radius) + '\n')
        # Add the peripheric material, alone on its line (geometric limit is
        # the pin pitch)
        out.write('mix' + str(materials[-1]) + '\n\n')
    # Writing the pin lattice
    out.write('lat 110 1 0.0 0.0 ' + str(npin) + ' ' + str(npin) + ' '
              + str(pinpitch) + '\n')
    for y in range(0,npin):
        pinline = [str(int(pin)) for pin in pinmap[:][y]]
        out.write(' '.join(pinline) + '\n')
    # Deal with water gap all around the assembly, composed of the same water
    # that is used for border cells
    cellName = cellNames[int(pinmap[0][1] - 1)]
    mat = fuel_geo[cellName]['MIX'][-1]
    # Python computation of npin*pinpitch/2 has only 8 precise digits. We must
    # truncate it. Otherwise the geometry becomes very slightly larger than it
    # should be, causing an 'undefined geometry' error in Serpent.
    out.write('\nsurf  1000  sqc  0.0 0.0 ' + str(truncate(npin*pinpitch/2,8))
              + '\n')
    out.write('surf  1001  sqc  0.0 0.0 ' + str(assemblypitch/2) + '\n')
    out.write('cell 110  0  fill 110   -1000' + '\n')
    out.write('cell 111  0  mix' + str(mat) + '       1000 -1001' + '\n')
    out.write('cell 112  0  outside     1001\n\n')
    #---
    #  Compositions
    #---
    compos = [] # List of compo objects
    library = loadlcm(filepath + '.compo')
    mix = library['ISOTOPESMIX']
    dens = library['ISOTOPESDENS']
    temp = library['ISOTOPESTEMP']
    name = library['ISOTOPERNAME'].split()
    for imix in np.unique(mix):
        compos.append(compo(filename,imix,temp,getMat(imix,mix,dens,name)))
    #---
    #  Write compositions
    #---
    for mix in mix_in_geom:
        for comp in compos:
            if mix == comp.getMix():
                comp.writeSerpent(out)
    out.close()
