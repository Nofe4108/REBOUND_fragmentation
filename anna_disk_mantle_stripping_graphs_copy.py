#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 22:15:13 2022

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

composition_input_file = "mantle_stripping_input/differentiated_mantle_stripping_input1.txt"
composition_output_file = "mantle_stripping_output/custom_mantle_stripping_output1.txt"
collision_report_file = "new_collision_reports/new_collision_report1.txt"
final_orbital_parameters_file = "final_orbital_parameters/final_orbital_parameters1.txt"
fig1_file = 'graphs/collision_type_bar_graph1.pdf'
fig2_file = 'graphs/final_core_fracs1.pdf'
fig3_file = 'graphs/initial_disk1.png'
fig4_file = 'graphs/final_planets1.png'

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

def calc_core_radius(mass, core_frac, core_density):
    core_mass = mass*core_frac
    core_radius = ((3*core_mass)/(4*core_density))**(1/3)
    return(core_radius)

def calc_radius(mass, core_frac, mantle_density, core_radius):
    mantle_frac = 1-core_frac
    mantle_mass = mass*mantle_frac
    radius = (core_radius**3+((3*mantle_mass)/(4*mantle_density)))**(1/3)
    return(radius)
      


###### MAIN FUNCTION #########
f = open(composition_output_file, 'r')
init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
compositions = organize_compositions(init_compositions) #organzies the data from the file
for line in compositions: #loops through each particle in the array
    if sum(line[2:4]) != 1.0: #makes sure composition fractions add up to one for each particle (uses second value in list and onwards)
        print ('ERROR: Realative abundances do not add up to 1.0')
        sys.exit(1)
                
    try: 
        init_hashes = [int(x[0]) for x in init_compositions] 
    except:
        init_hashes = [x[0].value for x in init_compositions]
    
file = open(collision_report_file, 'r')
blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list

collision_types = ['elastic bounce','merger','partial accretion','partial erosion','super catastrophic']
no_collisions = [0, 0, 0, 0, 0]
no_bounces = 0    
no_mergers = 0
no_accretions = 0
no_erosions = 0
no_catastrophics = 0

final_masses_em = []
final_core_fracs_em = []

for i in compositions:
    final_masses_em.append(i[1]*334672.021419) #converts to earth masses
    final_core_fracs_em.append(i[2])
for i in compositions:
    if .1 > i[1]*334672.021419 and i[1]*334672.021419  > .05:
        print(i[0])
for i in range(len(blocks)): #iterates through each value in blocks list - THIS IS A VERY BIG LOOP THAT CONTAINS ALL OF THE MASS TRANSFER DECISION MAKING
    block = blocks[i].split() #seperates each long string full of the collision data into its own list to be parsed through
    time = float(block[0]) #time of collision is first value
    collision_type = int(block[1]) #type of collision is second
    if collision_type == 0:
        no_collisions[0] += 1
    elif collision_type == 1:
        no_collisions[1] += 1
    elif collision_type == 2:
        no_collisions[2] += 1
    elif collision_type == 3:
        no_collisions[3] += 1
    else:
        no_collisions[4] += 1
        

fig1, ax1 = plt.subplots()
ax1.barh(collision_types, no_collisions)

plt.savefig(fig1_file, bbox_inches='tight', pad_inches=1.25)


#plt.show()


f = open(composition_input_file, 'r')
init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
original_compositions = organize_compositions(init_compositions) #organzies the data from the file
original_masses = [original_compositions[i][1]*334672.021419*1200 for i in range(len(original_compositions))]
original_core_fracs = [original_compositions[i][2] for i in range(len(original_compositions))]
original_a = [original_compositions[i][4] for i in range(len(original_compositions))]
original_e = [original_compositions[i][5] for i in range(len(original_compositions))]

f = open(final_orbital_parameters_file, 'r')
final_orbital_parameters = [line.split() for line in f.readlines()]
final_hashes = [int(final_orbital_parameters[i][0]) for i in range(len(final_orbital_parameters))]
final_masses = [float(final_orbital_parameters[i][1])*334672.021419*1200 for i in range(len(final_orbital_parameters))]
final_a = [float(final_orbital_parameters[i][2]) for i in range(len(final_orbital_parameters))]
final_e = [float(final_orbital_parameters[i][3]) for i in range(len(final_orbital_parameters))]
final_core_fracs = []
for hsh in final_hashes:
    for i in range(len(compositions)):
        if hsh == compositions[i][0]:
            final_core_fracs.append(compositions[i][2])

min_frac_init = min(original_core_fracs)
max_frac_init = max(original_core_fracs)
min_frac_final = min(final_core_fracs)
max_frac_final = max(final_core_fracs)

min_frac = 0
max_frac = 0
if min_frac_init <= min_frac_final:
    min_frac += min_frac_init
else:
    min_frac += min_frac_final
if max_frac_init >= max_frac_final:
    max_frac += max_frac_init
else:
    max_frac += max_frac_final

uniform_core_fracs = [.3 for i in range(len(original_core_fracs))]

fig2, ax2 = plt.subplots()
ax2.scatter(final_masses_em, final_core_fracs_em, color = 'black')
ax2.axhline(0.3, label = 'Initial Core Fraction', color = 'tab:blue', linestyle = '--', alpha=.7)
ax2.axvline(.093, label = 'Initial Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax2.axvline(.0093, label = 'Initial Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax2.set_xlabel('Mass ($M_{\u2295}$)',fontsize='large')
ax2.set_ylabel("Core Fraction", fontsize='large')
plt.xscale('log')
ax2.set_ylim([0, 1.1])
plt.yticks(np.arange(0, 1.1, .1))
ax2.minorticks_on()
plt.grid()
plt.legend()

plt.savefig(fig2_file, bbox_inches='tight', pad_inches=1.25)

fig3, ax3 = plt.subplots(figsize=(8,5))

color_map = plt.get_cmap('jet_r')
plot3 = ax3.scatter(original_a, original_e, marker='o', s=original_masses, c=original_core_fracs, cmap=color_map, vmin=min_frac, vmax=max_frac)
ax3.set_xlabel('Semi-major Axis (AU)',fontsize='large')
ax3.set_ylabel("Eccentricity", fontsize='large')
ax3.set_xlim(0.0, 4.5)
ax3.set_ylim(-0.01, 0.40)
cbar3 = fig3.colorbar(plot3, location='right', anchor=(0,0.5), pad=0.0)
cbar3.set_label(label='CMF', size='large')
cbar3.minorticks_on()
plt.grid()
#.0093*600
#.093*160

#Line2D([], [], color='black', marker='o', lw=0.0, label='0.01 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.01*1000))
 
ax3_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*1200))]                                                           

ax3_legend = plt.legend(handles=ax3_legend_elements, loc = 'upper right', prop={"size": 10})

plt.savefig(fig3_file, bbox_inches='tight', pad_inches=0.25, dpi=250)

fig4, ax4 = plt.subplots(figsize=(8,5))

plot4 = ax4.scatter(final_a, final_e, marker='o', s=final_masses, c=final_core_fracs, cmap=color_map, vmin=min_frac, vmax=max_frac)
ax4.set_xlabel('Semi-major Axis (AU)',fontsize='large')
ax4.set_ylabel("Eccentricity", fontsize='large')
ax4.set_xlim(0.0, 4.5)
ax4.set_ylim(-0.01, 0.40)
cbar4 = fig4.colorbar(plot4, location='right', anchor=(0,0.5), pad=0.0)
cbar4.set_label(label='CMF', size='large')
cbar4.minorticks_on()
plt.grid()
ax4_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*1000)),
                          Line2D([], [], color='black', marker='o', lw=0.0, label='0.5 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(0.5*1000))]  

ax4_legend = plt.legend(handles=ax3_legend_elements, loc = 'upper right', prop={"size": 10})

plt.savefig(fig4_file, bbox_inches='tight', pad_inches=0.25, dpi=250)


initial_mass = sum(original_masses)/1200
final_mass = sum(final_masses)/1200
non_uni_total_core_mass = 0
uni_total_core_mass = 0
final_total_core_mass = 0
avg_initial_core_frac = sum(original_core_fracs)/len(original_core_fracs)
for i in range(len(original_masses)):
    non_uni_total_core_mass += (original_masses[i]*original_core_fracs[i])/1200
    uni_total_core_mass += (original_masses[i]*0.3)/1200
    

for i in range(len(final_masses)):
    final_total_core_mass += (final_masses[i]*final_core_fracs[i])/1200

print(len(original_masses))
print(initial_mass)
print(final_mass)
print(non_uni_total_core_mass)
print(uni_total_core_mass)
print(final_total_core_mass)
print(avg_initial_core_frac)





