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
f = open("/Users/nferich/GitHub/REBOUND_fragmentation/disk_simulations/disk_mantle_stripping_input.txt", 'r')
init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
compositions = organize_compositions(init_compositions) #organzies the data from the file
for line in compositions: #loops through each particle in the array
    if sum(line[2:]) != 1.0: #makes sure composition fractions add up to one for each particle (uses second value in list and onwards)
        print ('ERROR: Realative abundances do not add up to 1.0')
        sys.exit(1)
                
    try: 
        init_hashes = [int(x[0]) for x in init_compositions] 
    except:
        init_hashes = [x[0].value for x in init_compositions]
    
file = open("/Users/nferich/GitHub/REBOUND_fragmentation/disk_simulations/collision_report.txt", 'r')
blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list

collision_types = ['elastic bounce','merger','partial accretion','partial erosion','super catastrophic']
no_collisions = [0, 0, 0, 0, 0]
no_bounces = 0    
no_mergers = 0
no_accretions = 0
no_erosions = 0
no_catastrophics = 0

final_masses = []
final_core_fracs = []

for i in compositions:
    final_masses.append(i[1]*334672.021419) #converts to earth masses
    final_core_fracs.append(i[2])

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


fig2, ax2 = plt.subplots()
ax2.scatter(final_masses, final_core_fracs, color = 'black')
ax2.axhline(0.3, label = 'Initial Core Fraction', color = 'tab:blue', linestyle = '--', alpha=.7)
ax2.axvline((6.4171e23/1.9885e30)*334672.021419, label = 'Initial Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax2.axvline((7.342e22/1.9885e30)*334672.021419, label = 'Initial Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax2.set_xlabel('Mass ($M_{\u2295}$)',fontsize='large')
ax2.set_ylabel("Core Fraction", fontsize='large')
plt.grid()
plt.legend()


#plt.show()

plt.savefig('/Users/nferich/GitHub/REBOUND_fragmentation/disk_simulations/collision_type_bar_graph.pdf', bbox_inches='tight', pad_inches=1.25)
plt.savefig('/Users/nferich/GitHub/REBOUND_fragmentation/disk_simulations/final_core_fracs.pdf', bbox_inches='tight', pad_inches=1.25)

