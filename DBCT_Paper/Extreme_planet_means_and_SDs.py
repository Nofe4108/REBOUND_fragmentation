"""
Produces Stats for planets with extreme CMFs from non-uniform distributions
"""


import astropy.constants as constants
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as matplotlib
import numpy as np
import math
import time
import sys
import matplotlib.patches as mpatches
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

# Constants
distribution_list = ['2step', 'lin', 'exp'] #Used to create full pathways to output files
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth

# Lists that will be used to produce data
final_masses = []
final_core_fracs = []
compositions = []

# File Pathways
comp_output_file_pw = "DBCT_output/"
comp_output_file_name = "_DBCT_output"


file_range = np.arange(1,51,1) # Number of output files to extract data from

for distr in distribution_list:
    distribution_masses = []
    distribution_core_fracs = []
    for no in file_range: 
        comp_output_file = comp_output_file_pw + distr + comp_output_file_name + str(no) + ".txt"
        f = open(comp_output_file, 'r')
        init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
        comp = organize_compositions(init_compositions) #organzies the data from the file
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

final_cmf_sds = [np.std(cmfs) for cmfs in final_planet_cmfs]
final_cmf_medians = [np.median(cmfs) for cmfs in final_planet_cmfs]
final_cmf_means = [np.mean(cmfs) for cmfs in final_planet_cmfs]
initial_cmf_means = [0.294, 0.314, 0.326]

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

extreme_planets_beneath_cutoff = []

for i, distr in enumerate(final_planet_extreme_masses):
    planet_counter = 0
    for j, mass in enumerate(distr):
        if mass <= final_extreme_mass_means[i]:
            if final_planet_extreme_cmfs[i][j] < final_cmf_means[i]-final_cmf_sds[i] or final_planet_extreme_cmfs[i][j] > final_cmf_means[i]+final_cmf_sds[i]:
                planet_counter+=1
    extreme_planets_beneath_cutoff.append(planet_counter)

#print(final_cmf_means)
#print(final_cmf_sds)
#print(final_extreme_mass_means)
print(final_cmf_medians)
print(extreme_planets_beneath_cutoff)
print(len(final_planet_extreme_masses[2]))
print(len(final_planet_masses[0]))