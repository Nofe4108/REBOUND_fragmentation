"""
Produces Figure 14 from Ferich et al. (IN PREP)
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
distribution_list = ['3step', 'lin', 'exp'] #Used to create full pathways to output files

# Lists that will be used to produce the plot
final_masses = []
final_cmfs = []

# File Pathways
DBCT_file_pw = "DBCT_output/"
DBCT_file_name = "_DBCT_output"
CT_file_pw = "CT_C&S_2022_output/"
CT_file_name = "_CT_output"

file_range = np.arange(1,51,1) #Number of output files to extract data from

"""Goes through and extracts mass and CMF data that was created with the Differentiated
Body Composition Tracker and the Composition Tracker"""

for distr in distribution_list:
    DBCT_final_masses = []
    DBCT_final_cmfs = []
    CT_final_masses = []
    CT_final_cmfs = []
    for no in file_range: 
        DBCT_output_file = DBCT_file_pw + distr + DBCT_file_name + str(no) + ".txt"
        CT_output_file = CT_file_pw + distr + CT_file_name + str(no) + ".txt"
        DBCT_f = open(DBCT_output_file, 'r')
        CT_f = open(CT_output_file, 'r')
        DBCT_init_compositions = [line.split() for line in DBCT_f.readlines()] # Reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
        CT_init_compositions = [line.split() for line in CT_f.readlines()]
        DBCT_comp = organize_compositions(DBCT_init_compositions) # Organzies the data from the file
        CT_comp = organize_compositions(CT_init_compositions)
        for line in DBCT_comp:
            if line[2] > 1.0 or line[2] < 0.0: # Makes sure CMF isn't negative or greater than one
                print ('ERROR: CMF from DBCT does not a realistic value - file number:', str(no))
        for line in CT_comp:
            if line[2] > 1.0 or line[2] < 0.0:
                print ('ERROR: CMF from CT does not a realistic value - file number:', str(no))
        for i in range(len(DBCT_comp)):
            DBCT_final_masses.append(DBCT_comp[i][1]*mass_conversion)
            DBCT_final_cmfs.append(DBCT_comp[i][2])
            CT_final_masses.append(CT_comp[i][1]*mass_conversion)
            CT_final_cmfs.append(CT_comp[i][2])
    final_masses.append(DBCT_final_masses)
    final_masses.append(CT_final_masses)
    final_cmfs.append(DBCT_final_cmfs)
    final_cmfs.append(CT_final_cmfs)
    
    
# Figure 14
fig, ([ax0, ax2, ax4], [ax1, ax3, ax5]) = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(13,10))
fig.subplots_adjust(wspace=0.03, hspace=0.03)
axes = [ax0, ax1, ax2, ax3, ax4, ax5] # List will be used to help plot all subplots
initial_average_cmfs = [0.294, 0.294, 0.314, 0.314, 0.326, 0.326] # Initial Average CMFS for each distriubtion

for axis in axes:
    axis_index = axes.index(axis)
    axis.grid(alpha=0.7)
    axis.scatter(final_masses[axis_index], final_cmfs[axis_index], s=3.0, color='black', zorder=2)
    axis.axhline(initial_average_cmfs[axis_index], label = 'Initial Average CMF', color = 'firebrick', linestyle = '--', alpha=.7, zorder=1) 
    axis.axvline(.093, label = 'Embryo Mass', color = 'royalblue', linestyle = ':', alpha=.7, zorder=1)
    axis.axvline(.0093, label = 'Planetesimal Mass', color = 'goldenrod', linestyle = '-.', alpha=.7, zorder=1)

ax0.legend(framealpha=1.0, fontsize=8.0)
ax0.set_title('a)', fontsize='x-large')
ax2.set_title('b)', fontsize='x-large')
ax4.set_title('c)', fontsize='x-large')
ax0.set_ylabel("CMF (DBCT)", fontsize='x-large', labelpad=10.0)
ax1.set_xlabel('Mass ($M_{\u2295}$)', fontsize='x-large', labelpad=10.0)
ax3.set_xlabel('Mass ($M_{\u2295}$)', fontsize='x-large', labelpad=10.0)
ax5.set_xlabel('Mass ($M_{\u2295}$)', fontsize='x-large', labelpad=10.0)
ax1.set_ylabel("CMF (CT)", fontsize='x-large', labelpad=10.0)
plt.xscale('log')
plt.ylim(-0.01, 1.01)
plt.yticks(np.arange(0.0, 1.1, 0.1))
plt.savefig("graphs/DBCT_and_CT_Final_CMFs.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/DBCT_and_CT_Final_CMFs.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/DBCT_and_CT_Final_CMFs.png', bbox_inches='tight', dpi=300)