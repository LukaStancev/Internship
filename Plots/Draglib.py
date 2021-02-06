#
#  Plotting Draglib content
#  Usage : python3 Draglib.py
#  Author : V. Salino (IRSN), 02/2021
#

# Imports
# faulthandler for segmentation faults. Left here due to a heisenbug in lcm
import faulthandler; faulthandler.enable()
import lcm
import numpy as np
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# We may plot quantiles or energy-correlated random samples
graphtype = 'quantiles'
#graphtype = 'samples'

for iso in glob.glob('../../Njoy/TENDL/*'):
    # Retrieve isotope name as given in directories' names
    iso = iso.rsplit('/', 1)[1]
    print('iso = ' + iso)
    #---
    #  Initialize plot
    #---
    fig, axs = plt.subplots(nrows = 2, sharex='all',
                            gridspec_kw={'hspace': 0})
    #---
    #  Create link toward Draglib file
    #---
    command = ('ln -s ../../Njoy/TENDL/' + iso + '/draglib' + iso + ' _draglib'
               + iso + '.txt')
    os.system(command)
    # A few plotting parameters
    if iso.startswith('U23'):
        minsig = 0.1
        loclegend = 'upper right'
        ncol = 2
    else:
        minsig = 0.01
        loclegend = 'lower left'
        ncol = 3
    # Initialize the flag indicating presence/absence (resp. True/False) of at
    # least one Autolib data
    AutolibFlag = False
    #---
    #  Load data (most time spent here)
    #---
    draglib = lcm.new('LCM_INP', 'draglib' + iso + '.txt')
    os.system('rm -f _draglib' + iso + '.txt')
    draglib._impx = 0
    #---
    #  Establish the list of available random samples, each being contained in
    #  a directory
    #---
    dirnames = [dirnam for dirnam in draglib.keys() if iso in dirnam]
    dirnames.sort()
    #---
    #  Establish the set of available reactions. Threshold reactions may be
    #  present in some random samples, but absent in others (if its cross
    #  sections are too low). The list established here is on a "minimum"
    #  basis, i.e. present in every random samples.
    #---
    reactions = None
    for dirname in dirnames:
        isotopedir = draglib[dirname]
        SubTemp = isotopedir['SUBTMP0002']
        if not reactions:
            reactions = set(SubTemp.keys())
        else:
            reactions = reactions.intersection(SubTemp.keys())
    #---
    #  Establish the list of interesting reactions among available ones, for
    #  that isotope. The notable reactions set aside are: NG, NELAS, NFTOT,
    #  N3N, N4N, NNP, N2A
    #---
    wished_reactions = {'NTOT0', 'NUSIGF', 'NINEL', 'N2N', 'NA', 'NP', 'ND'}
    if iso.startswith('U23'):
        wished_reactions = wished_reactions.difference({'NA', 'NP', 'ND'})
    elif iso.startswith('Zr9'):
        wished_reactions = wished_reactions.difference({'NA'})
    reactions = list(reactions.intersection(wished_reactions))
    reactions.sort()
    for reaction in reactions:
        print('reaction = ' + reaction)
        #---
        #  Recover the color so that all the curves of a given reaction are of
        #  the same color
        #---
        color = next(axs[0]._get_lines.prop_cycler)['color']
        #---
        #  Retrieve data that shall be plotted
        #---
        # Energy limits (e.g. 99g, 172g, 281g, 295g, 315g, 361g, ...)
        Energies = draglib['ENERGY']
        relstdEnergies = Energies
        # List that shall contain every random samples
        XSs = []
        for dirname in dirnames:
            irand = int(dirname[-3:])
            isotopedir = draglib[dirname]
            # Use 2nd temperature
            SubTemp = isotopedir['SUBTMP0002']
            Temperature = isotopedir['TEMPERATURE']
            Temperature = str(int(Temperature[1]))
            # Print list of possible reactions
            if irand == 0:
                print(SubTemp.keys())
            # Retrieve cross section
            XS = SubTemp[reaction]
            # Combine Autolib data and energy limits (where available) with XS
            # and energy limits outside of that (where Autolib is unavailable)
            if ('BIN-' + reaction) in SubTemp.keys():
                AutolibFlag = True
                # Sampling-insensitive data (identical in every random samples)
                if irand == 0:
                    AutolibBins = isotopedir['BIN-NFS']
                    BinEnergies = isotopedir['BIN-ENERGY']
                    Autolib_beg = np.min(np.nonzero(AutolibBins))
                    Autolib_end = np.max(np.nonzero(AutolibBins)) + 1
                    Ener_bef, _, Ener_aft = (np.split(Energies,
                                             [Autolib_beg, Autolib_end + 1]))
                    Energies = np.concatenate((Ener_bef, BinEnergies,
                                               Ener_aft))
                XS_bef, _, XS_aft = np.split(XS, [Autolib_beg, Autolib_end])
                Autolib = SubTemp['BIN-' + reaction]
                XS = np.concatenate((XS_bef, Autolib, XS_aft))
            else:
                # Minor (threshold) reactions may have less than G+1 cross
                # sections. The missing cross sections shall be considered as
                # equal to zero. Here, we strip the most thermal energy limits
                # in order to avoid plotting cross sections equal to zero.
                if (len(XS) + 1) < len(Energies):
                    Energies = Energies[:(len(XS) + 1)]
                    relstdEnergies = relstdEnergies[:(len(XS) + 1)]
            XSs.append(XS)
            if graphtype == 'samples':
                if irand == 0:
                    axs[0].step(Energies, np.append(XS, XS[-1]), where='post',
                                linewidth=0.05, color=color, label=reaction)
                elif irand < 20:
                    axs[0].step(Energies, np.append(XS, XS[-1]), where='post',
                                linewidth=0.05, color=color)
            # Deleting progressively to avoid segfaults in lcm module
            del SubTemp
            del isotopedir
        #---
        #  For some minor (threshold) reactions, we may come up with an
        #  irregular 2D array, because some samples may count more (non-zero)
        #  cross sections than other samples. We shall homogenize this by
        #  reducing to the smallest length.
        #---
        minlength = min(map(len, XSs))
        XSs = [XS[:minlength] for XS in XSs]
        #---
        #  Transform XSs into 2D NumPy arrays and then compute mean and
        #  relative standard deviation (in %) over every samples
        #---
        XSs = np.array(XSs)
        XS = np.mean(XSs, axis=0)
        relstdXS = np.std(XSs, axis=0)/XS*100
        #---
        #  Keep relative standard deviation only when the cross section is
        #  larger than minsig b (for threshold reactions)
        #---
        if XS[0] > minsig:
            largeXS = np.where(XS > minsig)[0]
            limit = len(np.intersect1d(largeXS, np.arange(0, len(largeXS))))
            relstdXS = relstdXS[:limit]
            relstdEnergies = Energies[:(len(relstdXS) + 1)]
        #---
        #  Add relative standard deviation plot
        #---
        axs[1].step(relstdEnergies, np.append(relstdXS, relstdXS[-1]),
                    where='post', linewidth=0.1, color=color)
        #---
        #  Add median and a few other quantiles
        #---
        if graphtype == 'quantiles':
            for q in [0.05, 0.25, 0.5, 0.75, 0.95]:
                qXS = np.quantile(XSs, q, axis=0)
                if q == 0.5:
                    axs[0].step(Energies, np.append(qXS, qXS[-1]),
                                where='post', linewidth=0.1, color=color,
                                label=reaction)
                else:
                    axs[0].step(Energies, np.append(qXS, qXS[-1]),
                                where='post', linewidth=0.1, color=color,
                                alpha=1/(abs(q-0.5)*3+1))
        #---
        #  Compute some correlations' informations
        #---
        if(reaction == 'NTOT0'):
            min_rij = np.min(np.corrcoef(XSs))
    del draglib
    #---
    #  Graph glitter
    #---
    # Set a title and axis labels
    axs[0].set_title('Statistics on ' + iso + ' cross section at '
                     + Temperature
                     + 'K versus energy, from TENDL Draglib', wrap=True)
    axs[1].set_xlabel('Incident neutron energy [eV]')
    if graphtype == 'samples':
        axs[0].set_ylabel('Cross section [b]')
    elif graphtype == 'quantiles':
        axs[0].set_ylabel('Cross section [b]\n(5%, 25%, 50%, 75%,\n '
                          + '95% quantiles)')
    axs[1].set_ylabel('Relative standard deviation')
    # Add legends
    leg = axs[0].legend(loc = loclegend, ncol = ncol)
    # Set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(0.6)
    # Cross sections are only readable in log-log graphs
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    # Add energy limits for Autolib domain
    if AutolibFlag:
        axs[0].axvline(BinEnergies[0], linestyle='--', color='black',
                       linewidth=0.2)
        axs[0].axvline(BinEnergies[-1], linestyle='--', color='black',
                       linewidth=0.2)
        axs[1].axvline(BinEnergies[0], linestyle='--', color='black',
                       linewidth=0.2)
        axs[1].axvline(BinEnergies[-1], linestyle='--', color='black',
                       linewidth=0.2)
    # Restrict the x-axis to data of interest. In particular, remove widely on
    # the left side : the most thermal group is enormous (~3 decades) and of
    # little interest.
    axs[0].set_xlim([Energies[-1]/10+Energies[-2]*9/10, Energies[0]])
    # Show only cross sections strictly greater than minsig b
    ylim = axs[0].get_ylim()
    axs[0].set_ylim(bottom = max(ylim[0], minsig*1.000001))
    # Start relative standard deviation at zero, to better see very low values
    # (at resonant peaks)
    axs[1].set_ylim(bottom = 0)
    # Relative standard deviation plot is capped at 90%
    ylimstd = min((axs[1].get_ylim())[1], 90)
    axs[1].set_ylim(top = ylimstd)
    # Show all powers of 10 in Energies xaxis, with 10 minor ticks
    locmaj = ticker.LogLocator(base = 10.0, numticks = 24)
    locmin = ticker.LogLocator(base = 10.0, numticks = 24,
                               subs = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
    axs[0].xaxis.set_major_locator(locmaj)
    axs[0].xaxis.set_minor_locator(locmin)
    axs[0].xaxis.set_minor_formatter(ticker.NullFormatter())
    # Set xaxis ticks on top instead of bottom
    axs[0].xaxis.set_ticks_position('top')
    # Add a light grid
    axs[0].grid(which='both', alpha=0.2, linewidth=0.1)
    axs[1].grid(which='both', alpha=0.2, linewidth=0.1)
    if ylimstd < 50:
        # Tick relative standard deviation every 5% (major) and 2.5% (minor)
        tickstd = 5.0
    else:
        # Tick relative standard deviation every 10% (major) and 5% (minor)
        tickstd = 10.0
    axs[1].yaxis.set_major_locator(ticker.MultipleLocator(tickstd))
    axs[1].yaxis.set_minor_locator(ticker.MultipleLocator(tickstd/2))
    # Set relative standard deviation y-ticks as ('5%',) '10%', ... except the
    # first one, which is a plain '0' (no percent sign)
    ytickspercent = ['0', '0']
    for i in np.arange(1, 19):
        ytickspercent.append(str(i*int(tickstd)) + '%')
    axs[1].set_yticklabels(ytickspercent)
    # Show minimum correlation coefficient in a corner of the figure
    axs[0].text(0.05, 0.9, r'$\mathrm{min}(\mathcal{R}_{ij})=$'
                + str(round(min_rij,5)), fontsize=8,
                transform=axs[0].transAxes)
    # Tight layout prevents a label from exceeding the figure frame
    fig.tight_layout()
    #---
    #  Save plot as pdf (vectorized)
    #---
    fig.savefig('XS_' + iso + '_' + graphtype + '.pdf')
    #---
    #  Clean-up for next plot
    #---
    plt.close('all')
    del fig
    del axs
os.system('gs -q -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=XS_' + graphtype
          + '.pdf -dBATCH XS_*_' + graphtype + '.pdf')
print("Plotting completed")
