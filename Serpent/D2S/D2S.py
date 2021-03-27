#  Convert library-type and geometry-type DRAGON objects to Serpent input file.
#  It is limited to PWR geometries symmetrical by eighth.
#  Usage   : from D2S import D2S
#            InfiniteLattice(geofilename)
#            or
#            FullCore(geofilename)
#  Authors : V. Salino, B. Dechenaux Sensei
#  Date    : 11/2020

#---
# Imports
#---
import math
import numpy as np
import os
import decimal
from collections import OrderedDict
# faulthandler for segmentation faults. Left here due to a heisenbug in lcm.
import faulthandler; faulthandler.enable()
import lcm

xsdata = '../PyNjoy2016/output/Universal.xsdata'
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
    def writeSerpent(self, output, mixsuffix):
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
                moder = 'moder lwtr' + str(self.m_mix) + mixsuffix + ' 1001\n'
        # Write the header line for the entire material...
        output += 'mat mix' + str(self.m_mix) + mixsuffix + ' sum\n'
        output += moder
        output += 'tmp ' + str(self.m_temp) + ' % Kelvin\n'
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
            with open('../' + xsdata) as xsdatafile:
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
                    aceFile = xsdata_subset[i-1][0]
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
            output += "%s %.8E\n" %(aceFile,self.m_compo[iso])
        # Write the 'thermr' card, if needed
        if tslFiles:
            output += ('therm lwtr' + str(self.m_mix) + mixsuffix + ' '
                       + str(self.m_temp) + ' ' + tslFiles + '\n')
        output += '\n'
        return output
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
def sss2Assembly(filepath, filename, pinprefix, mixsuffix, univassmb):
    # Universe number attributed to this specific assembly
    #  = 0 for infinite lattice models
    # != 0 for full core models
    univassmb = str(univassmb)
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
    #---
    #  Write assembly geometry
    #---
    output = ''
    # Writing the few different type of pincells
    for cellID in range(0,len(cellNames)):
        cellName = cellNames[cellID]
        radii = fuel_geo[cellName]['RADIUS'][1:]
        materials = fuel_geo[cellName]['MIX']
        # Writes cellID together with its name as a commentary
        output += 'pin ' + pinprefix + str(cellID+1) + ' %' + cellName + '\n'
        for radius,material in zip(radii, materials):
            output += ('mix' + str(material) + mixsuffix + ' ' + str(radius)
                       + '\n')
        # Add the peripheric material, alone on its line (its geometric limit
        # is the pin pitch)
        output += 'mix' + str(materials[-1]) + mixsuffix + '\n\n'
    # Writing the pin lattice
    output += ('lat ' + str(900 + int(univassmb)) + ' 1 0.0 0.0 ' + str(npin)
               + ' ' + str(npin) + ' ' + str(pinpitch) + '\n')
    for y in range(0,npin):
        pinline = [pinprefix + str(int(pin)) for pin in pinmap[:][y]]
        output += ' '.join(pinline) + '\n'
    # Deal with water gap all around the assembly, composed of the same water
    # that is used for border cells
    cellName = cellNames[int(pinmap[0][1] - 1)]
    mat = fuel_geo[cellName]['MIX'][-1]
    # Python computation of npin*pinpitch/2 has only 8 precise digits. We must
    # truncate it. Otherwise the geometry becomes very slightly larger than it
    # should be, causing an 'undefined geometry' error in Serpent.
    output += ('\nsurf  ' + str(1000 + int(univassmb)) + '  sqc  0.0 0.0 '
               + str(truncate(npin*pinpitch/2,8)) + '\n')
    output += ('surf  ' + str(1100 + int(univassmb)) + ' sqc  0.0 0.0 '
               + str(assemblypitch/2) + '\n')
    output += ('cell ' + str(100 + int(univassmb)) + '  ' + univassmb
               + '  fill ' + str(900 + int(univassmb)) + '   -'
               + str(1000 + int(univassmb)) + '\n')
    output += ('cell ' + str(200 + int(univassmb)) + '  ' + univassmb
               + '  mix' + str(mat) + mixsuffix + '       '
               + str(1000 + int(univassmb)) + ' -'
               + str(1100 + int(univassmb)) + '\n')
    if univassmb == 0:
        output += ('cell ' + str(300 + int(univassmb)) + '  ' + univassmb
                   + '  outside     ' + str(1100 + int(univassmb)) + '\n\n')
    #---
    #  Assess the compositions used in this geometry
    #---
    mix_in_geom = set()
    for cellID in range(0,len(cellNames)):
        cellName = cellNames[cellID]
        materials = fuel_geo[cellName]['MIX']
        mix_in_geom.update(materials.tolist())
    mix_in_geom = sorted(list(mix_in_geom))
    #---
    #  Retrieve the compositions
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
                output = comp.writeSerpent(output, mixsuffix)
    return output, assemblypitch
# trimRefl function removes, in a deterministic 1D reflector geometry:
# * the fuel zone included in reflector models,
# * the fictious limits, i.e. when neighbouring meshes have identical mixes.
def trimRefl(mix, meshx):
    trimmed_mix = []
    trimmed_meshx = []
    if mix[-1] == 2:
        raise Exception('MIX 2 is placed in last position and should not')
    for i in range(0, len(mix)):
        # Recognize MIX 2 as equivalent to the one that follows it
        # It cannot overflow, thanks to the raise Exception above
        if mix[i] == 2:
            mix[i]=mix[i+1]
        # In a reflector model, MIX 1 means fuel. Do not retain it.
        if mix[i] != 1:
            # Embark the first one, whatever it is
            if len(trimmed_mix) == 0:
                trimmed_mix.append(mix[i])
                trimmed_meshx.append(meshx[i])
            else:
                # Embark only if it is different from the last one embarked
                # (same as previous indicates a fictious limit)
                if mix[i] != trimmed_mix[-1]:
                    trimmed_mix.append(mix[i])
                    trimmed_meshx.append(meshx[i])
    # Embark last limit
    trimmed_meshx.append(meshx[-1])
    return trimmed_mix,trimmed_meshx
# Main function for 2D infinice lattices
def InfiniteLattice(filepath):
    global compo
    filename = filepath.split('/')[-1]
    #---
    #  Prepare file writing for Serpent
    #---
    sss2 = open(filename + '.sss2', 'w')
    #---
    #  Write general information
    #---
    sss2.write('% Serpent 2 dataset produced automatically with D2S.\n')
    sss2.write('set title ' + '"Tihange ' + filename + '"\n')
    sss2.write('set acelib "' + xsdata + '"\n')
    sss2.write('set bc periodic\n')
    sss2.write('set pop 6000 500 20\n')
    sss2.write('set ures 1\n')
    sss2.write('set edepmode 2\n')
    sss2.write('set gcu -1\n')
    sss2.write('set cmm 0\n')
    sss2.write('set shbuf 0 0\n')
    sss2.write('set repro 0\n')
    sss2.write('plot 3 2500 2500\n')
    sss2.write('mesh 3 2500 2500\n\n')
    #
    output = sss2Assembly(filepath, filename, '', '', 0)
    sss2.write(output)
    sss2.close()
# Main function for 3D full cores
def FullCore(filepath, Assemblies, CB):
    fullcore_geo = loadlcm(filepath)
    #---
    #  Write general information
    #---
    sss2 = open('TihangeFullCore' + str(CB) + 'ppm.sss2', 'w')
    sss2.write('% Serpent 2 dataset produced automatically with D2S.\n')
    sss2.write('set title ' + '"Tihange full core, ' + str(CB) + ' ppm"\n')
    sss2.write('set acelib "' + xsdata + '"\n')
    sss2.write('set bc black\n')
    sss2.write('set pop 500000 1000 200 1.0\n')
    sss2.write('set ures 1\n')
    sss2.write('set edepmode 2\n')
    sss2.write('set gcu -1\n')
    sss2.write('set cmm 0\n')
    sss2.write('set shbuf 0 0\n')
    sss2.write('set repro 0\n')
    sss2.write('\n')
    #---
    #  Perform a few consistency checks
    #---
    meshx = fullcore_geo['MESHX']
    meshy = fullcore_geo['MESHY']
    if len(meshx) != len(meshy):
        raise Exception('Non-square cores are not supported.')
    sizex = meshx[1:] - meshx[:-1]
    sizey = meshy[1:] - meshy[:-1]
    if not np.allclose(sizex, sizey):
        raise Exception('All meshes are not identical on the X-Y plane.')
    #---
    #  Retrieve loading pattern
    #---
    meshz = fullcore_geo['MESHZ']
    sizez = meshz[1:] - meshz[:-1]
    mixs = fullcore_geo['MIX']
    # Get it into 3D numpy array
    mixs = mixs.reshape(len(sizez), len(sizey), len(sizex))
    # Remove reflectors (=1, 2 or 3)
    mixs[mixs <= 3] = 0
    # Replace duplicate assemblies with their siblings
    UniqueAssemblies = {}
    Duplicates = {}
    for key, value in Assemblies.items():
        if value not in UniqueAssemblies.values():
            UniqueAssemblies[key] = value
        else:
            Duplicates[key] = value
            # If it is a duplicate, find its sibling assembly
            newkey = list(UniqueAssemblies.keys())[list(UniqueAssemblies.values()).index(value)]
            # Replace the duplicate value with its sibling
            mixs = np.where(mixs == key, newkey, mixs)
    # Remove empty and duplicate layers and remove in-between fictious limits
    # in meshz
    newmixs = []
    newsize = []
    for i in range(0,len(sizez)):
        # Do not take into account empty layers (composed only of zeroes)
        if not ((len(np.unique(mixs[i])) == 1) and (mixs[i, 0, 0] == 0)):
            identical = False
            if i > 0:
                identical = not np.array_equal(mixs[i], mixs[i-1])
            if (i == 0) or identical:
                newmixs.append(mixs[i, 1:-1, 1:-1])
                newsize.append(sizez[i])
            else:
                # If a layer is not counted, we shall add its thickness to the
                # thickness of the previous (identical) layer
                newsize[-1] = newsize[-1] + sizez[i]
    mixs = np.array(newmixs)
    sizez = np.array(newsize)
    # Now that the empty layers are removed, we can put back a non-zero
    # universe all around the core.
    mixs[mixs == 0] = 1
    # Write uniquely every assembly in that core, except the surroundings ('1')
    for iassembly in (set(np.unique(mixs)) - {1}):
        assemblyfilepath = ('/'.join(filepath.split('/')[:-1]) + '/UOX'
                            + Assemblies[iassembly] + '_' + str(CB) + 'ppm')
        assemblyfilename = filepath.split('/')[-1]
        sss2.write("% Assembly " + Assemblies[iassembly] + '\n\n')
        #
        output, assemblypitch = sss2Assembly(assemblyfilepath,
                                             assemblyfilename,
                                             f"{iassembly:02d}",
                                             '_' + Assemblies[iassembly],
                                             iassembly)
        sss2.write(output)
    sss2.write('% Full core\n\n')
    sss2.write('surf 2 inf\n')
    # Void around the core (will be replaced by reflector geometry)
    sss2.write('cell 3 1 void -2\n')
    # The assemblies' pitch from DONJON is often not the most representative :
    # instead, we shall use the DRAGON assembly pitch, as it includes adequate
    # thermal expansion
    sss2.write('surf 4 sqc 0.0 0.0 ' + str(assemblypitch*(len(mixs[0,0,:]))/2)
               + '\n')
    # Layers separating different axial fuel zones
    for i in range(0, len(sizez) - 1):
        sss2.write('surf ' + str(5 + i) + ' pz ' + str(sum(sizez[:i+1]))
                   + '\n')
    sss2.write('\n')
    # Write loading pattern
    for i in range(0, len(sizez)):
        sss2.write('lat ' + str(30 + i) + ' 1 0.0 0.0 ' + str(len(sizex) - 2)
                   + ' ' + str(len(sizey) - 2) + ' ' + str(assemblypitch)
                   + '\n')
        for y in range(0, len(mixs[i, :, 0])):
            assemblyline = [str(assembly) for assembly in mixs[i, :, y]]
            sss2.write(' '.join(assemblyline) + '\n')
        # Fill this layer with loading pattern layer
        axial_layers = ''
        if i == 0:
            axial_layers += '1400' # Bottom of the core
        else:
            axial_layers += str(4 + i)
        if i == (len(sizez) - 1):
            axial_layers += ' -1600' # Top of the core
        else:
            axial_layers += ' -' + str(5 + i)
        sss2.write('cell ' + str(70 + i) + ' 0 fill ' + str(30 + i) + ' -1200 '
                   + axial_layers + '\n\n')
    #---
    #  Radial reflector geometry
    #---
    sss2.write('% Radial reflector\n')
    radial_geo = loadlcm('../geo_compo/ReflRadial_' + str(CB) + 'ppm.geo')
    mix, meshx = trimRefl(radial_geo['MIX'], radial_geo['MESHX'])
    # Steel baffle from the lower right quarter of the core
    gcross = []
    middle = int((len(mixs[0, :, :]) - 1)/2)
    for assemblyline in mixs[0, middle+1:, middle+1:]:
        nassembly = np.count_nonzero(assemblyline > 1)
        gcross.append(np.count_nonzero(assemblyline > 1))
    # Remove duplicates
    gcross = list(OrderedDict.fromkeys(gcross))
    # Compute the baffle distances
    gcross = (np.array(gcross)+0.5)*assemblypitch
    sss2.write('\nsurf 1200 gcross 0.0 0.0 '
               + ' '.join([str(i) for i in gcross]))
    gcross = gcross + (meshx[1] - meshx[0])
    sss2.write('\nsurf 1201 gcross 0.0 0.0 '
               + ' '.join([str(i) for i in gcross]))
    sss2.write('\n')
    # After the baffle, a series of cylinders is used (thermal shield, etc)
    for i in range(2, len(meshx)):
        sss2.write('surf ' + str(1200 + i) + ' cyl 0.0 0.0 '
                   + str(meshx[i]) + '\n')
    for i in range(0, len(mix)):
        sss2.write('\ncell ' + str(1300 + i) + ' 0 mix' + str(mix[i])
                   + 'radial ' + str(1200 + i) + ' -' + str(1201 + i)
                   + ' 1400 -1600')
    # The outside is beyond that
    sss2.write('\ncell ' + str(1300 + len(mix)) + ' 0 outside '
               + str(1200 + len(mix)) + ' 1400 -1600')
    del radial_geo
    #---
    #  Radial reflector composition
    #---
    # Assess the compositions used in this geometry
    mix_in_geom = sorted(list(set(mix)))
    # Retrieve the compositions
    compos = [] # List of compo objects
    radial_compo = loadlcm('../geo_compo/ReflRadial_' + str(CB) + 'ppm.compo')
    mix =  radial_compo['ISOTOPESMIX']
    dens = radial_compo['ISOTOPESDENS']
    temp = radial_compo['ISOTOPESTEMP']
    name = radial_compo['ISOTOPERNAME'].split()
    for imix in np.unique(mix):
        compos.append(compo('',imix,temp,getMat(imix,mix,dens,name)))
    output = ''
    # Write compositions
    for mix in mix_in_geom:
        for comp in compos:
            if mix == comp.getMix():
                output = comp.writeSerpent(output, 'radial')
    sss2.write('\n\n' + output)
    del radial_compo
    #---
    #  Bottom reflector geometry
    #---
    sss2.write('% Bottom reflector\n')
    bottom_geo = loadlcm('../geo_compo/ReflBottom_' + str(CB) + 'ppm.geo')
    mix, meshb = trimRefl(bottom_geo['MIX'], bottom_geo['MESHX'])
    for i in range(0, len(meshb)):
        sss2.write('\nsurf ' + str(1400 + i) + ' pz ' + str(-meshb[i]))
    sss2.write('\n')
    for i in range(0, len(mix)):
        sss2.write('\ncell ' + str(1500 +i) + ' 0 mix' + str(mix[i])
                   + 'bottom -' + str(1400 + i) + ' ' + str(1401 +i))
    sss2.write('\ncell ' + str(1500 + len(mix)) + ' 0 outside -'
               + str(1400 + len(mix)))
    #---
    #  Bottom reflector composition
    #---
    # Assess the compositions used in this geometry
    mix_in_geom = sorted(list(set(mix)))
    # Retrieve the compositions
    compos = [] # List of compo objects
    bottom_compo = loadlcm('../geo_compo/ReflBottom_' + str(CB) + 'ppm.compo')
    mix =  bottom_compo['ISOTOPESMIX']
    dens = bottom_compo['ISOTOPESDENS']
    temp = bottom_compo['ISOTOPESTEMP']
    name = bottom_compo['ISOTOPERNAME'].split()
    for imix in np.unique(mix):
        compos.append(compo('',imix,temp,getMat(imix,mix,dens,name)))
    output = ''
    # Write compositions
    for mix in mix_in_geom:
        for comp in compos:
            if mix == comp.getMix():
                output = comp.writeSerpent(output, 'bottom')
    sss2.write('\n\n' + output)
    del bottom_compo
    #---
    #  Top reflector geometry
    #---
    sss2.write('% Top reflector\n')
    top_geo = loadlcm('../geo_compo/ReflTop_' + str(CB) + 'ppm.geo')
    mix, mesht = trimRefl(top_geo['MIX'], top_geo['MESHX'])
    for i in range(0, len(mesht)):
        sss2.write('\nsurf ' + str(1600 + i)
                   + ' pz ' + str(mesht[i] + sum(sizez)))
    sss2.write('\n')
    for i in range(0, len(mix)):
        sss2.write('\ncell ' + str(1700 +i) + ' 0 mix' + str(mix[i]) + 'top '
                   + str(1600 + i) + ' -' + str(1601 +i))
    sss2.write('\ncell ' + str(1700 + len(mix)) + ' 0 outside '
               + str(1600 + len(mix)))
    #---
    #  Top reflector composition
    #---
    # Assess the compositions used in this geometry
    mix_in_geom = sorted(list(set(mix)))
    # Retrieve the compositions
    compos = [] # List of compo objects
    top_compo = loadlcm('../geo_compo/ReflTop_' + str(CB) + 'ppm.compo')
    mix =  top_compo['ISOTOPESMIX']
    dens = top_compo['ISOTOPESDENS']
    temp = top_compo['ISOTOPESTEMP']
    name = top_compo['ISOTOPERNAME'].split()
    for imix in np.unique(mix):
        compos.append(compo('',imix,temp,getMat(imix,mix,dens,name)))
    output = ''
    # Write compositions
    for mix in mix_in_geom:
        for comp in compos:
            if mix == comp.getMix():
                output = comp.writeSerpent(output, 'top')
    sss2.write('\n\n' + output)
    del top_compo
    #---
    #  To account for possible energy deposition in the surroundings, tally on
    #  a grid larger than the core:
    #  - radially, with 2 more assemblies in each direction,
    #  - axially, with 3 bins : below the core, the active core itself, and
    #    above the core. Each of these 3 axial bins has the same size as the
    #    active core.
    #---
    xyn=len(mixs[0,0,:]) + 2*2
    xymax=str(assemblypitch*xyn/2)
    xyn=str(xyn)
    det = ('\n'
           + 'dx -' + xymax + ' ' + xymax +  ' ' + xyn + '\n'
           + 'dy -' + xymax + ' ' + xymax +  ' ' + xyn + '\n'
           + 'dz -' + str(sum(sizez)) + ' ' + str(2*sum(sizez)) + ' 3\n')
    zpixel = str(int(2500*meshx[-1]/(meshb[-1] + sum(sizez) + mesht[-1])))
    #---
    #  Write plots and tallies options
    #---
    sss2.write('% Plots and tallies options\n')
    sss2.write('plot 3 5000 5000 100.0\n')
    sss2.write('plot 2 ' + zpixel + ' 2500\n')
    sss2.write('mesh 3 5000 5000\n')
    sss2.write('det -4' + det + 'dr -4 void\n')
    sss2.write('det -6' + det + 'dr -6 void\n')
    sss2.write('det -7' + det + 'dr -7 void\n')
    sss2.write('det -8' + det + 'dr -8 void\n')
    sss2.write('det -80' + det + 'dr -80 void\n')
    sss2.write('set cpd 1 1 0.0 ' + str(sum(sizez)) + '\n')
    sss2.close()
