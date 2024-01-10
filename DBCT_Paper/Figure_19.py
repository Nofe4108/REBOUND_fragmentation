"""
Produces Figure 19 from Ferich et al. (IN PREP)
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth

#Lists that will be used to produce the plot
final_hashes = []
final_masses = []
final_cmfs = []

# File pathways
comp_output_file_pw = "DBCT_output/2step_DBCT_output"

file_range = np.arange(1,51,1) #Number of output files to extract data from

for no in file_range:
    composition_output_file = comp_output_file_pw + str(no) + ".txt"
    f = open(composition_output_file, 'r')
    init_compositions = [line.split() for line in f.readlines()] # Reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) # Organzies the data from the file
    for obj in compositions:
        final_hashes.append(obj[0])
        final_masses.append(obj[1]*mass_conversion)
        final_cmfs.append(obj[2])

# Figure 19
plt.rcParams['axes.axisbelow'] = True # Makes sure grid is behind points
fig1, ax1 = plt.subplots(figsize=(6,5))
ax1.scatter(final_masses, final_cmfs, s=5.0, color = 'black')
ax1.axhline(0.304, label = 'Initial CMF', color = 'tab:blue', linestyle = '--', alpha=.7)
ax1.axvline(.093, label = 'Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax1.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax1.set_xlabel('Mass ($M_{\u2295}$)',fontsize='large')
ax1.set_ylabel("CMF", fontsize='large')
plt.xscale('log')
ax1.set_ylim([-0.01, 1.1])
plt.yticks(np.arange(0, 1.1, .1))
ax1.minorticks_on()
ax1.grid(alpha=0.7)
plt.legend()
plt.savefig("graphs/fig19.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig19.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig19.png', bbox_inches='tight', dpi=300)