#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 17:12:34 2023

@author: nferich
"""

from astropy import units as u
import astropy.constants as constants
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import math
import time
import sys

core_density = 7874.0 #kg m^-3 - density of iron
mantle_density = 3330.0 #kg m^3 -  approximate density of Earth's upper mantle

composition_input_file1 = "mantle_stripping_input/nu_mantle_stripping_input1.txt"
composition_input_file2 = "mantle_stripping_input/lin_mantle_stripping_input1.txt"
composition_input_file3 = "mantle_stripping_input/exp_mantle_stripping_input1.txt"
fig_file = 'graphs/initial_disks.pdf'


######## COMPOSITION DATA ORGANIZING FUNCTION ###########
#Function takes the data from the composition input file and organizes it in a list
#The created list has seperate values that represent the particle's hash, mass, core frac, and mantle frac
def organize_compositions(init_compositions):
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

f = open(composition_input_file1, 'r')
init_compositions1 = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
original_compositions1 = organize_compositions(init_compositions1) #organzies the data from the file
original_masses1 = [original_compositions1[i][1]*334672.021419*1200 for i in range(len(original_compositions1))]
original_core_fracs1 = [original_compositions1[i][2] for i in range(len(original_compositions1))]
original_a1 = [original_compositions1[i][3] for i in range(len(original_compositions1))]
original_e1 = [original_compositions1[i][4] for i in range(len(original_compositions1))]

f = open(composition_input_file2, 'r')
init_compositions2 = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
original_compositions2 = organize_compositions(init_compositions2) #organzies the data from the file
original_masses2 = [original_compositions2[i][1]*334672.021419*1200 for i in range(len(original_compositions2))]
original_core_fracs2 = [original_compositions2[i][2] for i in range(len(original_compositions2))]
original_a2 = [original_compositions2[i][3] for i in range(len(original_compositions2))]
original_e2 = [original_compositions2[i][4] for i in range(len(original_compositions2))]

f = open(composition_input_file3, 'r')
init_compositions3 = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
original_compositions3 = organize_compositions(init_compositions3) #organzies the data from the file
original_masses3 = [original_compositions3[i][1]*334672.021419*1200 for i in range(len(original_compositions3))]
original_core_fracs3 = [original_compositions3[i][2] for i in range(len(original_compositions3))]
original_a3 = [original_compositions3[i][3] for i in range(len(original_compositions3))]
original_e3 = [original_compositions3[i][4] for i in range(len(original_compositions3))]

min_frac_init = min(original_core_fracs3)
max_frac_init = max(original_core_fracs3)


#fig3, ax3 = plt.subplots(figsize=(7,5))

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, gridspec_kw={'width_ratios': [10, 10, 10.5]}, figsize=(10,5))
fig.subplots_adjust(wspace=0)
ax1.minorticks_on()
#ax2.grid()
#ax3.grid()
ax1.set_title('a)', fontsize='x-large')
ax2.set_title('b)', fontsize='x-large')
ax3.set_title('c)', fontsize='x-large')
fig.supxlabel('Semi-major Axis (AU)', fontsize='x-large')
ax1.set_ylabel("Eccentricity", fontsize='x-large', labelpad=10.0)
plt.xlim(0.0, 4.5)
plt.ylim(0.000, 0.0101)

color_map = plt.get_cmap('jet')

p1 = ax1.scatter(original_a1, original_e1, marker='o', s=original_masses1, c=original_core_fracs1, cmap=color_map, vmin=min_frac_init, vmax=max_frac_init+.01)
p2 = ax2.scatter(original_a2, original_e2, marker='o', s=original_masses2, c=original_core_fracs2, cmap=color_map, vmin=min_frac_init, vmax=max_frac_init+.01)
p3 = ax3.scatter(original_a3, original_e3, marker='o', s=original_masses3, c=original_core_fracs3, cmap=color_map, vmin=min_frac_init, vmax=max_frac_init+.01)


cbar = fig.colorbar(p3, location='right', anchor=(0,0), pad=0.0, ticks=[.1, .2, .3, .4, .5, .6, .7])
cbar.set_label(label='CMF', size='x-large', labelpad=10.0)
cbar.minorticks_on()

plt.savefig(fig_file, bbox_inches='tight', pad_inches=0.01)
