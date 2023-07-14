#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 13:45:45 2023

@author: nferich
"""

from astropy import units as u
import astropy.constants as constants
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as matplotlib
import numpy as np
import math
import time
import sys
import matplotlib.patches as mpatches


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

###### MAIN FUNCTION #########


final_masses = []
final_core_fracs = []
compositions = []

comp_output_file_pw = "mantle_stripping_output/"
comp_output_file_name = "_mantle_stripping_output"


file_range = np.arange(1,51,1)
mass_conversion = 334672.021419
distribution_list = ['nu', 'lin', 'exp']

############## START OF LOOP #########################
for distr in distribution_list:
    distribution_masses = []
    distribution_core_fracs = []
    for no in file_range: 
        comp_output_file = comp_output_file_pw + distr + comp_output_file_name + str(no) + ".txt"
        f = open(comp_output_file, 'r')
        init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
        comp = organize_compositions(init_compositions) #organzies the data from the file
        for line in comp: #loops through each particle in the array
            if line[2] > 1.0 or line[2] < 0.0: #makes sure CMF isn't negative or greater than one
                print ('ERROR: CMF does not a realistic value')
                sys.exit(1)             
                try: 
                    init_hashes = [int(x[0]) for x in init_compositions] 
                except:
                    init_hashes = [x[0].value for x in init_compositions]
        for obj in comp:
            distribution_masses.append(obj[1]*mass_conversion)
            distribution_core_fracs.append(obj[2])
    final_masses.append(distribution_masses)
    final_core_fracs.append(distribution_core_fracs)
 
final_planet_masses = []
final_planet_cmfs = []


for i in range(len(final_masses)):
    bp_masses = []
    bp_cmfs = []
    for j in range(len(final_masses[i])):
            if final_masses[i][j] >= 0.093:
                bp_masses.append(final_masses[i][j])
                bp_cmfs.append(final_core_fracs[i][j])
    final_planet_masses.append(bp_masses)
    final_planet_cmfs.append(bp_cmfs)

final_cmf_sds = [np.std(final_planet_cmfs[i]) for i in range(len(final_planet_masses))]
final_cmf_medians = [np.median(final_planet_cmfs[i]) for i in range(len(final_planet_masses))]
final_cmf_means = [np.mean(final_planet_cmfs[i]) for i in range(len(final_planet_masses))]
initial_cmf_means = [0.294, 0.314, 0.326]
print(final_cmf_means)

final_planet_extreme_masses = []
final_planet_extreme_cmfs = []
for i in range(len(final_planet_masses)):
    masses = []
    cmfs = []
    for j in range(len(final_planet_masses[i])):
        if final_planet_cmfs[i][j] < final_cmf_means[i]-final_cmf_sds[i] or final_planet_cmfs[i][j] > final_cmf_means[i]+final_cmf_sds[i]:
            masses.append(final_planet_masses[i][j])
            cmfs.append(final_planet_cmfs[i][j])
    final_planet_extreme_masses.append(masses)
    final_planet_extreme_cmfs.append(cmfs)
            
final_extreme_mass_means = [np.mean(final_planet_extreme_masses[i]) for i in range(len(final_planet_extreme_masses))]
final_extreme_cmf_means = [np.mean(final_planet_extreme_cmfs[i]) for i in range(len(final_planet_extreme_cmfs))]
print(final_cmf_sds)
print(final_extreme_mass_means)
print(final_cmf_medians)
