"""
Produces the data used in Table 2 in 
Ferich et al. (IN PREP)
"""

import numpy as np
from astropy import units as u
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
distribution_list = ['3step', 'lin', 'exp']

# Lists that will be used to produce numbers for the table - Contains obects that are considered planets
final_masses = []
final_cmfs = []
compositions = []

# File pathways
composition_output_file_pw = "DBCT_output/"
composition_output_file_name = "_DBCT_output"

file_range = np.arange(1,51,1) #Number of output files to extract data from

############## START OF LOOP #########################
for distr in distribution_list:
    bp_masses = [] # bp stands for big planet
    bp_cmfs = []
    for no in file_range: 
        composition_output_file = composition_output_file_pw + distr + composition_output_file_name + str(no) + ".txt"
        f = open(composition_output_file, 'r')
        init_compositions = [line.split() for line in f.readlines()] # Reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
        f.close()
        comp = organize_compositions(init_compositions) # Organzies the data from the file
        for obj in comp:
            if obj[1]*mass_conversion >= 0.093:
                bp_masses.append(obj[1]*mass_conversion)
                bp_cmfs.append(obj[2])
    final_masses.append(bp_masses)
    final_cmfs.append(bp_cmfs)
 

final_cmf_sds = [np.std(distr_cmfs) for distr_cmfs in final_cmfs]
final_cmf_medians = [np.median(distr_cmfs) for distr_cmfs in final_cmfs]
final_cmf_means = [np.mean(distr_cmfs) for distr_cmfs in final_cmfs]
initial_cmf_means = [0.294, 0.314, 0.326]
final_extreme_masses = []
final_extreme_cmfs = []

for i in range(len(final_masses)):
    extreme_masses = []
    extreme_cmfs = []
    for j in range(len(final_masses[i])):
        if final_cmfs[i][j] < final_cmf_means[i]-final_cmf_sds[i] or final_cmfs[i][j] > final_cmf_means[i]+final_cmf_sds[i]:
            extreme_masses.append(final_masses[i][j])
            extreme_cmfs.append(final_cmfs[i][j])
    final_extreme_masses.append(extreme_masses)
    final_extreme_cmfs.append(extreme_cmfs)
            
final_extreme_mass_means = [np.mean(distr_extreme_masses[i]) for distr_extreme_masses in final_extreme_masses]
final_extreme_cmf_means = [np.mean(distr_extreme_cmfs[i]) for distr_extreme_cmfs in final_extreme_cmfs]

print(final_cmf_sds)
print(final_extreme_mass_means)
print(final_cmf_medians)

"""planets_above_cutoff = [[] for i in range(len(final_masses))]
for i, distr in distribution_list:
    for j, obj_mass in enumerate(bp_masses[i]):
        if obj_mass > final"""
        
   
    