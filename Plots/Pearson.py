#
#  Compute Pearson correlation coefficient (squared) between:
#  * power of the central assembly, and
#  * hydrogen (H1) cross sections
#  Usage: python3 Pearson.py
#  Author: V. Salino (IRSN), 02/2021
#
# -*- coding: utf-8 -*-

# Imports
import lcm
import numpy as np
from scipy.stats import pearsonr
import os
import glob
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from BasicFunctions import *
from matplotlib.ticker import ScalarFormatter
matplotlib.rcParams['text.usetex'] = False


plt.rcParams.update(tex_fonts())
lang = 'fr'  # fr/en

# Retrieve H1 cross sections, through (very basic) ENDF parsing

from GetXSE import *
os.system('./ParseH1H2O.sh')
GetISO()
choice = input("Select the isotope: ")
Energies, XSs, reactions, iso, isotopes, keys_to_find = GetXSE(choice)
print('Energies', Energies)
print('size Energies', np.shape(Energies))
for key in keys_to_find:
    print('samples', key)

# Initialize a gaussian vector (used for legends)
gaussian = np.random.randn(300)
# Declare control rod insertions and their full names
controlrods = ['CD', 'D', 'ARO']
#controlrods = ['D']
text = {}
text['ARO'] = 'all rods out'
text['D'] = 'D rod bank inserted'
text['CD'] = 'C and D rod banks inserted'
# Create a dictionary to store Pearson Correlation Coefficient and its pvalue
# between the central assembly and the H1_H2O cross section
PCC = {}
pvalue = {}
XSs_new = {}
energies = {}
reactions_to_skip = {}
for controlrod in controlrods:
    print('controlrod = ' + controlrod)
    # Create links toward every available power distribution
    os.system('ln -s ../Drakkar/Output_TIH_TMC/_Power' + controlrod
              + '_*.ascii .')
    # List isotopes we've been (potentially) randomly sampling
    firstfile = glob.glob('_Power' + controlrod + '_*.ascii')[0]
    ResultFile = lcm.new('LCM_INP', firstfile[1:])
    isotopes = ResultFile['NamIso'].split()

    del ResultFile
    # Initialize dictionary (isotope-wise) of lists (sampling-wise) that will
    # contain the retrieved data
    powers = {}
    for iso in isotopes:
        powers[iso] = []
    isamples = []
    for file in list(glob.glob('_Power' + controlrod + '_*.ascii')):
        # Loading data while removing the beginning-of-filename underscore
        ResultFile = lcm.new('LCM_INP', file[1:])
        os.system("rm " + file)
        # Compute mean to normalize power to 1.0 per average assembly
        power = ResultFile['POWER-CHAN']
        mean = np.sum(power) / len(power)
        power = power / mean
        # Determine which isotope was randomly sampled to store the
        # obtained power distribution in the right place
        for iso in isotopes:
            if (int(ResultFile[iso]) != -33) and (iso == choice):
                powers[iso].append(power)
                # Also recover the sample number, to order them and compute
                # correlations
                isamples.append(int(ResultFile[choice]))
        del ResultFile
    # Update 'isotopes' list so that only those isotopes that were indeed
    # randomly sampled can be the subject of loops to come
    isotopes_used = []
    for iso in isotopes:
        if len(powers[iso]) != 0:
            isotopes_used.append(iso)
    isotopes = isotopes_used
    # Order them to facilitate the calculation of correlations
    a1 = []
    for iso in isotopes:
        print('current :' + iso)
        tmp = []
        
        for (a, b) in sorted(zip(isamples, powers[iso])):            
            tmp.append(b)
            a1.append(a)
        powers[iso] = tmp
    print('isamples',a1)
    # Convert this list of 1D numpy array into a 2D numpy array:
    # [irand, iassembly]
    tmp = {}
    for iso in isotopes:
        tmp[iso] = np.array(powers[iso])
    powers = tmp
    # Update 'isotopes' list so that only those isotopes that were indeed
    # randomly sampled can be the subject of loops to come
    isotopes_used = []
    for iso in isotopes:
        if powers[iso].size != 0:
            isotopes_used.append(iso)
    isotopes = isotopes_used
    #Update XSs to coresponding power 
    all_index = np.arange(300)
         #calculate missing indexes
    missing_index = np.setdiff1d(all_index, a1)
         # Updateing XSs array
    XSs_new[controlrod] = {}
    for reaction in XSs:
        XSs_new[controlrod][reaction] = np.delete(XSs[reaction], missing_index, axis=0)
        
    # Layout of the assemblies in the core
    CoreLayout, FullCoreLayout = GetCoreLayout(len(powers[isotopes[0]][0, :]))
    #
    FullPowers2D = {}
    Powers2D = {}
    for iso in isotopes:
        # Deploy the 1D-indexed power maps into 2D-indexed power maps
        FullPowers2D[iso] = np.zeros((len(powers[iso][:, 0]),  # irand
                                      len(FullCoreLayout[:, 0]),  # Y-axis
                                      len(FullCoreLayout[0, :])  # X-axis
                                      ))
        i = 0
        for x in range(0, len(FullCoreLayout[0, :])):
            for y in range(0, len(FullCoreLayout[:, 0])):
                if FullCoreLayout[x, y] != 0:
                    FullPowers2D[iso][:, x, y] = powers[iso][:, i]
                    i = i + 1
        # Fold in quarter maps (averaging per rotation) for plotting purposes
        Powers2D[iso] = np.zeros((len(powers[iso][:, 0]),  # irand
                                  len(CoreLayout[:, 0]),  # Y-axis
                                  len(CoreLayout[0, :])  # X-axis
                                  ))
        for irand in range(0, len(FullPowers2D[iso][:, 0, 0])):
            Powers2D[iso][irand, :, :] = fold(FullPowers2D[iso][irand, :, :])

    # Compute Pearson's correlation coefficient (PCC) for the central assembly
    PCC[controlrod] = {}
    pvalue[controlrod] = {}
    energies[controlrod] = {}
    energies_temp = Energies[:-1]
    reactions_to_skip[controlrod] = []
    for reaction in XSs_new[controlrod]:
		# Create an energy array that corresponds to each reaction
        if reaction.startswith('N'): 
           if len(np.shape(XSs_new[controlrod][reaction])) == 1:
              reactions_to_skip[controlrod].append(reaction)
              continue
           elif np.shape(XSs_new[controlrod][reaction])[1] == len(energies_temp):
                energies[controlrod][reaction] = np.array(energies_temp)
           elif np.shape(XSs_new[controlrod][reaction])[1] < len(energies_temp):
                column_length = np.shape(XSs_new[controlrod][reaction])[1]
                energies[controlrod][reaction] = np.empty(column_length)
                energies[controlrod][reaction][:column_length] = energies_temp[:column_length]
                energies[controlrod][reaction][column_length:] = energies_temp[-1]
           else:
               print(f"Skipping reaction {reaction} in control rod {controlrod} due to IndexError in shape.")
               continue
                      
        print('reaction', reaction)
        if reaction in reactions_to_skip:
            print('reaction skiped', reaction)
            continue
			
        if reaction.startswith('N') and np.shape(XSs_new[controlrod][reaction]) == (len(Powers2D[choice][:, 7, 0]),len(energies[controlrod][reaction])):
            print('reaction included', reaction)
            PCC[controlrod][reaction] = np.zeros_like(energies[controlrod][reaction])
            pvalue[controlrod][reaction] = np.zeros_like(energies[controlrod][reaction])            
            for i in np.arange(0,len(energies[controlrod][reaction])):
                # Focus on the central assembly
                PCC[controlrod][reaction][i], pvalue[controlrod][reaction][i] = pearsonr(XSs_new[controlrod][reaction][:, i], Powers2D[choice][:, 7, 0])
        else:
            # If the shape is not as expected, skip this reaction
            continue
            PCC[controlrod][reaction] = np.full_like(energies[controlrod][reaction], np.nan)
            pvalue[controlrod][reaction] = np.full_like(energies[controlrod][reaction], np.nan)
for controlrod in controlrods:
    for reaction in reactions_to_skip[controlrod]:
        del XSs_new[controlrod][reaction]

# Correlations plots
reactions_to_plot = ['NELAS', 'NG', 'NTOT0', 'ND', 'NINEL', 'NP', 'NA', 'N2N', 'N2N', 'NUSIGF']

for controlrod in ['CD', 'D', 'ARO']:
#for controlrod in ['D']:
    fig, ax = plt.subplots()    
    for reaction in XSs_new[controlrod]:
        if reaction.startswith('N') and np.shape(XSs_new[controlrod][reaction]) == (len(Powers2D[choice][:, 7, 0]),len(energies[controlrod][reaction])):
            if reaction == 'NELAS':
                txt = '(n,el)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NG':
                txt = r'(n,$\gamma$)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NTOT0':
                txt = '(n,tot)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'ND':
                txt = '(n,d)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NINEL':
                txt = '(n,inel)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NP':
                txt = '(n,p)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NA':
                txt = '(n,$\\alpha$)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'N2N':
                txt = '(n,2n)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NUSIGF':
                txt = '(n,NUSIGF)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', linestyle='--',  label=txt)
            elif reaction == 'N2A':
                txt = '(n,2$\\alpha$)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NPP':
                txt = '(n,2p)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)   
            elif reaction == 'N2N':
                txt = '(n,2n)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'N3N':
                txt = '(n,3n)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)  
            elif reaction == 'N4N':
                txt = '(n,4n)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NFTOT':
                txt = '(n,ftot)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NJJS00':
                txt = '(n,NJJS00)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
            elif reaction == 'NNP':
                txt = '(n,np)'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)                           
            else:
                txt = 'void'
                ax.step(energies[controlrod][reaction], PCC[controlrod][reaction]**2, where='post', label=txt)
        else:
               continue
        
            

    # Graph glitter
    if lang == 'en':
        ax.set_title('Correlations between the uncertainty of the central assembly\npower and the elastic and capture cross sections of hydrogen,\n' + text[controlrod])
        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('Pearson correlation coefficient, squared')
    elif lang == 'fr':
        ax.set_xlabel('Énergie [eV]')
        ax.set_ylabel('Coefficient de corrélation de Pearson au carré')
    ax.set_xscale('log')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlim(left=min(energies_temp), right=max(energies_temp[:-1]))
    ax.set_ylim(bottom=0, top=1)
    

    #ax.grid(which='both', alpha=0.2, linewidth=0.1)
    

    # Save plot as pdf (vectorized)
    os.system('mkdir -p output_Pearson_Draglib_Power/' + choice)
    fig.savefig('output_Pearson_Draglib_Power/' + choice + '/' + controlrod + '.pdf', bbox_inches='tight')
    fig.savefig('output_Pearson_Draglib_Power/' + choice + '/' + controlrod + '.svg', bbox_inches='tight')
    plt.close(fig)

# Scatter plots
fig, ax1 = plt.subplots()

# Scatter plot for 'nelastic' reaction
if 'NELAS' in XSs and np.shape(XSs_new['ARO']['NELAS']) == (len(Powers2D[choice][:, 7, 0]),len(energies[controlrod][reaction])):
    txt1 = '(n,el)'
    ax1.scatter(XSs_new['ARO']['NELAS'][:, 191], Powers2D[choice][:, 7, 0], label=txt1, marker='x', color='blue')
    if lang == 'en':
        ax1.set_xlabel('Cross section (n,el)')
        ax1.set_ylabel('Power (n,el)')
    elif lang == 'fr':
        ax1.set_xlabel('Section efficace (n,el) à 1 eV ')
        ax1.set_ylabel("Puissance de l'assemblage central (n,el)")
    ax1.set_xscale('log')
    ax1.set_xlim(left=min(XSs['NELAS'][:, 191]), right=max(XSs['NELAS'][:, 191]))
    ax1.legend(loc='upper left')
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    ax1.xaxis.set_minor_formatter(ScalarFormatter())
    ax1.ticklabel_format(axis='x', style='plain')

    # Create a second y-axis for 'ngamma' reaction if it exists
    if 'NG' in XSs and np.shape(XSs['NG']) == (len(Powers2D[choice][:, 7, 0]),len(energies[controlrod][reaction])):
        ax2 = ax1.twinx()
        txt2 = r'(n,$\\gamma$)'
        ax2.scatter(XSs_new['ARO']['NELAS'][:, 191], Powers2D[choice][:, 7, 0], label=txt2, marker='x', color='red')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel("Puissance de l'assemblage central (n,$\gamma$)")
        ax2.legend(loc='upper right')

        # Create a second x-axis for 'ngamma' reaction
        ax3 = ax1.twiny()
        ax3.scatter(XSs_new['ARO']['NELAS'][:, 191], Powers2D[choice][:, 7, 0], marker='x', color='red')
        ax3.xaxis.tick_top()
        ax3.xaxis.set_label_position("top")
        ax3.set_xlabel('Section efficace (n,$\gamma$) à 1 eV ')

        # Adjust the position of the second x-axis labels
        ax3.xaxis.set_label_coords(0.5, 1.1)

        # Adjust the position of the second x-axis tick labels
        ax3.xaxis.set_tick_params(pad=10)

    # Set custom ticks and spacing for the y-axis
    ax1.set_yticks(ax1.get_yticks()[::2])
    ax1.yaxis.set_tick_params(pad=5)

    # Set custom ticks and spacing for the x-axis
    ax1.set_xticks(ax1.get_xticks()[::2])
    ax1.xaxis.set_tick_params(pad=5)

    # Save plot as pdf (vectorized)
    os.makedirs('output_Pearson_Draglib_Power/' + choice, exist_ok=True)
    #fig.savefig(f'output_Pearson_Draglib_Power/{choice}/Power_Xss_1eV_common_Luka.pdf', bbox_inches='tight')
    plt.close(fig)
else:
    plt.close(fig)
    
# Clean-up
plt.close('all')

