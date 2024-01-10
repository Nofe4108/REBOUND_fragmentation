# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 18:11:51 2023

@author: nferi

Produces Figure 9 and 10 from Ferich et al. (IN PREP). 
Also produces the data inputted into Table 1
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from astropy import units as u
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
ef_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:red'] # Colors that represent each expansion factor from sim: blue-3, orange-5, green-7, purple-10, red-15

############ DATA COLLECTION #################

#Lists that will be used to produce final plot - bp stands for Big Planet - ef stands for expansion factor
final_hashes = []
final_masses = []
final_cmfs = []
efs = []
bp_final_hashes = [[] for i in range(len(ef_colors))]
bp_final_masses = [[] for i in range(len(ef_colors))]
bp_final_cmfs = [[] for i in range(len(ef_colors))]

comp_output_file_pw = "DBCT_output/uni_DBCT_output"
file_range = np.arange(1,51,1) #Number of output files to extract data from

for no in file_range:
    composition_output_file = comp_output_file_pw + str(no) + ".txt"
    f = open(composition_output_file, 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) #organzies the data from the file
    for obj in compositions:
        final_hashes.append(obj[0])
        final_masses.append(obj[1]*mass_conversion)
        final_cmfs.append(obj[2])
        if no < 11:
            efs.append('tab:blue')
        elif 11 <= no < 21:
            efs.append('tab:orange')
        elif 21 <= no < 31:
            efs.append('tab:green')
        elif 31 <= no < 41:
            efs.append('tab:purple')
        elif 41 <= no < 51:
            efs.append('tab:red')


        
# Extracts masses and compositional data for planets produced by REBOUND simulations
for i in range(len(final_masses)):
    if final_masses[i] >= 0.093:
        for j in range(len(ef_colors)):
            if efs[i] == ef_colors[j]:
                bp_final_hashes[j].append(final_hashes[i])
                bp_final_masses[j].append(final_masses[i])
                bp_final_cmfs[j].append(final_cmfs[i])
            
########### PLOTTING ######################
    
plt.rcParams['axes.axisbelow'] = True # Makes sure grid is behind points

# Figure 9
fig1, ax1 = plt.subplots(figsize=(6,5))
n_bins = np.linspace(.25, .36, num=22, endpoint=False)
#n_bins = np.linspace(.05, .55, num=50, endpoint=False)
#n_bins = np.linspace(0.00, 0.55, num=55, endpoint=False)
hist_type = 'barstacked'
ax1.hist(bp_final_cmfs, bins=n_bins, histtype=hist_type, lw=3.0, color=ef_colors, alpha = 1.0)
ax1.set_xlim([0.24, 0.36])
ax1.set_xlabel('CMF',fontsize='large')
ax1.set_ylabel('Number of Planets (M>0.093 $M_{\u2295}$)', fontsize='large')
#plt.yticks(np.arange(0, 1.1, .1))
ax1.minorticks_on()
ax1.grid(alpha=0.7)
ax4_legend_elements = [mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:blue', facecolor='tab:blue', label='ef=3'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:orange', facecolor='tab:orange', label='ef=5'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:green', facecolor='tab:green', label='ef=7'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:purple', facecolor='tab:purple', label='ef=10'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:red', facecolor='tab:red', label='ef=15')]
legend = plt.legend(handles=ax4_legend_elements, loc = 'upper right', framealpha = 0.7)
plt.gca().add_artist(legend)
plt.savefig("graphs/fig9.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig9.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig9.png', bbox_inches='tight', dpi=300)
    
# Figure 10
n_bins = np.linspace(.25, .36, num=132, endpoint=False)
fig2, ax2 = plt.subplots(figsize=(6,5))
for i in range(len(bp_final_cmfs)):
    ax2.hist(bp_final_cmfs[i], bins=n_bins, density='True', histtype='step', lw=1.8, color=ef_colors[i], alpha = 0.7, cumulative='True')
ax2.set_xlabel('CMF',fontsize='large')
ax2.set_ylabel('Probability of Occurrence', fontsize='large')
ax2.grid(alpha=0.7)
ax2.minorticks_on()
ax2.set_xlim([0.25, 0.358])
ax2.set_ylim([-0.01, 1.05])
ax2_legend_elements = [mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:blue', facecolor='tab:blue', label='ef=3'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:orange', facecolor='tab:orange', label='ef=5'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:green', facecolor='tab:green', label='ef=7'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:purple', facecolor='tab:purple', label='ef=10'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:red', facecolor='tab:red', label='ef=15')]
legend = plt.legend(handles=ax2_legend_elements, loc = 'upper left', framealpha = .7)
plt.gca().add_artist(legend)
plt.savefig("graphs/fig10.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig10.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig10.png', bbox_inches='tight', dpi=300)
    
# Table 1
statistics = []
p_values = []
for i in(range(len(bp_final_cmfs))):
    stat = []
    p_val = []
    for j in range(len(bp_final_cmfs)):
        test_result = stats.kstest(bp_final_cmfs[i], bp_final_cmfs[j], alternative='two-sided')
        stat.append(test_result.statistic)
        p_val.append(test_result.pvalue)
    statistics.append(stat)
    p_values.append(p_val)
print(statistics)
print(p_values)

ef_values = ['3', '5', '7', '10', '15']
print("Expansion Factor KS Test Results:")
for i in range(len(p_values)):
    for j in range(len(p_values[i])):
        if i < j:
            print("P-value between", ef_values[i], "and", ef_values[j]+': ', str(p_values[i][j]))