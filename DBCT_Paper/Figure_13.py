# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 16:05:28 2023

@author: nferi

Produces Figure 13 from Ferich et al. (IN PREP)
Need a new version of the Organize Compositions 
functions to extract the data
"""

import matplotlib as matplotlib
import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

def organize_compositions(init_compositions):
    """Data organizing function
    
    Takes in raw string data extracted from mantle stripping input file
    Creates a list filled with properties of each particle
    The properties include hash, mass, and CMF
    
    Parameters:
    init_compositions (list) -- raw particle data from input file

    Returns:
    compositions (list) -- nested list with properly formatted particle data
    """
    compositions = []
    for line in init_compositions: #iterates through each line from file
        particle_data = [] #This will hold data for individual particle - will contain its hash, mass, and the fractions for each layer
        for i in range(len(line)): #goes through the index for each element in a line
            if i == 0: #First value in each line is the hash so that'll be an integer
                particle_data.append(int(line[i]))
            else: #The other values will be floats
                particle_data.append(float(line[i]))
        compositions.append(particle_data) #add all the data about the particle to the compositions list
    return(compositions)

composition_input_file1 = "DBCT_input/3step_DBCT_input.txt"
composition_input_file2 = "DBCT_input/lin_DBCT_input.txt"
composition_input_file3 = "DBCT_input/exp_DBCT_input.txt"

mass_conversion = 334672.021419 # Msun to Mearth

f = open(composition_input_file1, 'r')
init_compositions1 = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
original_compositions1 = organize_compositions(init_compositions1) #organzies the data from the file
original_masses1 = [original_compositions1[i][1]*334672.021419*1200 for i in range(len(original_compositions1))]
original_core_fracs1 = [original_compositions1[i][2] for i in range(len(original_compositions1))]
original_a1 = [original_compositions1[i][3] for i in range(len(original_compositions1))]
original_e1 = [original_compositions1[i][4] for i in range(len(original_compositions1))]
f.close()

f = open(composition_input_file2, 'r')
init_compositions2 = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
original_compositions2 = organize_compositions(init_compositions2) #organzies the data from the file
original_masses2 = [original_compositions2[i][1]*334672.021419*1200 for i in range(len(original_compositions2))]
original_core_fracs2 = [original_compositions2[i][2] for i in range(len(original_compositions2))]
original_a2 = [original_compositions2[i][3] for i in range(len(original_compositions2))]
original_e2 = [original_compositions2[i][4] for i in range(len(original_compositions2))]
f.close()

f = open(composition_input_file3, 'r')
init_compositions3 = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
original_compositions3 = organize_compositions(init_compositions3) #organzies the data from the file
original_masses3 = [original_compositions3[i][1]*334672.021419*1200 for i in range(len(original_compositions3))]
original_core_fracs3 = [original_compositions3[i][2] for i in range(len(original_compositions3))]
original_a3 = [original_compositions3[i][3] for i in range(len(original_compositions3))]
original_e3 = [original_compositions3[i][4] for i in range(len(original_compositions3))]
f.close()

min_frac_init = min(original_core_fracs3)
max_frac_init = max(original_core_fracs3)

# Figure 13
color_map = cmr.sapphire_r
fig1, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, gridspec_kw={'width_ratios': [10, 10, 10.5]}, figsize=(10,5))
fig1.subplots_adjust(wspace=0)
p1 = ax1.scatter(original_a1, original_e1, marker='o', s=original_masses1, c=original_core_fracs1, cmap=color_map, vmin=min_frac_init, vmax=max_frac_init+.01)
p2 = ax2.scatter(original_a2, original_e2, marker='o', s=original_masses2, c=original_core_fracs2, cmap=color_map, vmin=min_frac_init, vmax=max_frac_init+.01)
p3 = ax3.scatter(original_a3, original_e3, marker='o', s=original_masses3, c=original_core_fracs3, cmap=color_map, vmin=min_frac_init, vmax=max_frac_init+.01)
ax1.minorticks_on()
ax1.set_title('a)', fontsize='x-large')
ax2.set_title('b)', fontsize='x-large')
ax3.set_title('c)', fontsize='x-large')
fig1.supxlabel('Semi-major Axis (AU)', fontsize='x-large')
ax1.set_ylabel("Eccentricity", fontsize='x-large', labelpad=10.0)
plt.xlim(0.0, 4.5)
plt.ylim(0.000, 0.0101)
cbar = fig1.colorbar(p3, location='right', anchor=(0,0), pad=0.0, ticks=[.1, .2, .3, .4, .5, .6, .7])
cbar.set_label(label='CMF', size='x-large', labelpad=10.0)
cbar.minorticks_on()
plt.savefig("graphs/fig13.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig13.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig13.png', bbox_inches='tight', dpi=300)



