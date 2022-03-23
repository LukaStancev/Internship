import glob
import numpy as np

def LoadH1H2O():
    renergies = {}
    rXSs = {}
    for reaction in ['nelastic', 'ngamma']:
        samples = []
        isamples = []
        for thisfile in glob.glob('strip/*_' + reaction):
            filename = thisfile.split('/')[-1]
            isample = filename.split('-')[-1].split('_')[0]
            isamples.append(int(isample))
            with open(thisfile) as f:
                sample = f.readlines()
                # Remove whitespace characters like `\n` at the end of each
                # line
                sample = [i.strip() for i in sample]
                # Split energy and cross section value
                sample = [i.split() for i in sample]
                samples.append(sample)
        # Order the samples in order to compute correlations
        tmp = []
        for (a,b) in sorted(zip(isamples, samples)):
            tmp.append(b)
        samples = tmp
        # Verify that all the energy grids have the same length
        for sample in samples:
            if len(sample) != len(samples[0]):
                raise Exception('One sample appears to have a different '
                                + 'length')
        nbsamples = len(samples)
        nbepoints = len(samples[0])
        energies = []
        XS = []
        XSs = []
        for j in np.arange(0, nbepoints):
            for i in np.arange(0, nbsamples):
                # Verify that only one energy grid exists among the samples.
                # The comparison is performed on strings (from the ENDF file),
                # which is better than on floats because we can strictly test
                # their equality.
                if samples[i][j][0] != samples[0][j][0]:
                    raise Exception('One sample appears to have a different '
                                    + 'energy grid')
                # Retrieve the cross sections
                XS.append(float(samples[i][j][1]))
            XSs.append(XS)
            XS = []
            # Retrieve the energy grid
            energies.append(float(samples[i][j][0]))
        # Transform XSs and energies float lists into NumPy arrays
        renergies[reaction] = np.array(energies)
        rXSs[reaction] = np.array(XSs)
    return renergies, rXSs
