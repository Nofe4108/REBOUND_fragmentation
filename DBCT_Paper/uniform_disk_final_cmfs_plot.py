"""
Produces Figure ? from Ferich et al. (IN PREP)
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions
import cmasher as cmr
import random

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
#ef_colors = ['#E69F00', '#4BA1D2', '#11775B', '#F0E442', '#9C2202'] # Colors that represent each expansion factor from sim: orange-3, light_blue-5, green-7, yellow/dark_blue-10, red-15 #F0E442 for yellow #002175 for dark blue #4BA1D2 for light blue 
#ef_colors = ['#648FFF','#785EF0','#DC267F','#FE6100','#FFB000']
ef_colors = ['firebrick', 'darkcyan', 'goldenrod', 'royalblue', 'darkmagenta']        
ef_markers = ['o','s','d','^','*']
ef_marker_sizes = [7.0, 5.0, 7.0, 7.0, 10.0]

#Lists that will be used to produce the plot
final_hashes = [[] for i in range(len(ef_colors))]
final_masses = [[] for i in range(len(ef_colors))]
final_cmfs = [[] for i in range(len(ef_colors))]


# File Pathways
comp_output_file_pw = "DBCT_output/uni_DBCT_output"

file_range = np.arange(1,51,1) #Number of output files to extract data from

for no in file_range:
    composition_output_file = comp_output_file_pw + str(no) + ".txt"
    f = open(composition_output_file, 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) #organzies the data from the file
    for obj in compositions:
        if no < 11:
            final_hashes[0].append(obj[0])
            final_masses[0].append(obj[1]*mass_conversion)
            final_cmfs[0].append(obj[2])
        elif 11 <= no < 21:
            final_hashes[1].append(obj[0])
            final_masses[1].append(obj[1]*mass_conversion)
            final_cmfs[1].append(obj[2])
        elif 21 <= no < 31:
            final_hashes[2].append(obj[0])
            final_masses[2].append(obj[1]*mass_conversion)
            final_cmfs[2].append(obj[2])
        elif 31 <= no < 41:
            final_hashes[3].append(obj[0])
            final_masses[3].append(obj[1]*mass_conversion)
            final_cmfs[3].append(obj[2])
        elif 41 <= no < 51:
            final_hashes[4].append(obj[0])
            final_masses[4].append(obj[1]*mass_conversion)
            final_cmfs[4].append(obj[2])


combined_data_points = []
for i in range(len(ef_colors)):
    for j in range(len(final_masses[i])):
        data_point = [final_masses[i][j], final_hashes[i][j], final_cmfs[i][j], ef_colors[i], ef_markers[i], ef_marker_sizes[i]]
        combined_data_points.append(data_point)

random.Random(1).shuffle(combined_data_points)
                
# Figure 4
plt.rcParams['axes.axisbelow'] = True # Makes sure grid is behind points
fig1, ax1 = plt.subplots(figsize=(6,5))
#for i, masses in enumerate(final_masses):
   # ax1.scatter(masses, final_cmfs[i], s=7.0, c=ef_colors[i], marker=ef_markers[i], vmin=0.0, vmax=1.0, zorder=2, alpha=0.8)
for i, data_point in enumerate(combined_data_points):
    ax1.scatter(data_point[0], data_point[2], s=data_point[5], c=data_point[3], marker=data_point[4], zorder=2, alpha=0.8)
    
ax1.axhline(0.3, label = 'Initial CMF', color = 'black', linestyle = '--', alpha=0.5, zorder=1)
ax1.axvline(.093, label = 'Embryo Mass', color = 'black', linestyle = ':', alpha=0.5, zorder=1)
ax1.axvline(.0093, label = 'Planetesimal Mass', color = 'black', linestyle = '-.', alpha=0.5, zorder=1)
ax1.set_xlabel('Mass ($M_{\u2295}$)',fontsize='x-large')
ax1.set_ylabel("CMF", fontsize='x-large')
plt.xscale('log')
ax1.set_ylim([-0.01, 0.8])
plt.yticks(np.arange(0, 0.8, .1))
ax1.minorticks_on()
ax1.grid(alpha=0.7)
"""ax1_legend_elements = [mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor= ef_colors[0], facecolor=ef_colors[0], label='ef=3'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor= ef_colors[1], facecolor=ef_colors[1], label='ef=5'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor= ef_colors[2], facecolor=ef_colors[2], label='ef=7'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor= ef_colors[3], facecolor=ef_colors[3], label='ef=10'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor= ef_colors[4], facecolor=ef_colors[4], label='ef=15')]"""
ax1_legend_elements = [Line2D([], [], lw=0.0, color=ef_colors[0], marker=ef_markers[0], markersize=7.0, label='ef=3'),
                       Line2D([], [], lw=0.0, color=ef_colors[1], marker=ef_markers[1], markersize=7.0, label='ef=5'),
                       Line2D([], [], lw=0.0, color=ef_colors[2], marker=ef_markers[2], markersize=7.0, label='ef=7'),
                       Line2D([], [], lw=0.0, color=ef_colors[3], marker=ef_markers[3], markersize=7.0, label='ef=10'),
                       Line2D([], [], lw=0.0, color=ef_colors[4], marker=ef_markers[4], markersize=7.0, label='ef=15')] 
legend1 = plt.legend(handles=ax1_legend_elements, loc = 'upper right', framealpha = .7)
legend2 = plt.legend(loc='upper center')
plt.gca().add_artist(legend1)
plt.savefig("/Users/nferi/Downloads/uni_final_cmfs.png", bbox_inches='tight', dpi=500)