"""
Produces Figure 20 from Ferich et al. (IN PREP)
"""

import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
from astropy import units as u
from matplotlib.lines import Line2D
import random
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
distribution_list = ['3step', 'lin', 'exp'] #Used to create full pathways to output files
ef_markers = ['o','s','d','^','*']
ef_marker_sizes = [15.0, 9.0, 15.0, 15.0, 30.0]

# Lists that will be used for plotting
final_masses = []
final_cmfs = []
final_sa = []

# File pathways
composition_output_file_pw = "DBCT_output/"
composition_output_file_name = "_DBCT_output"
final_orbital_parameters_file_pw = "final_orbital_parameters/final_orbital_parameters"

file_range = np.arange(1,51,1)

for distr in distribution_list:
    distr_final_masses = [[] for i in range(len(ef_markers))]
    distr_final_cmfs = [[] for i in range(len(ef_markers))]
    distr_final_sa = [[] for i in range(len(ef_markers))]
    for no in file_range:
        composition_output_file = composition_output_file_pw + distr + composition_output_file_name + str(no) + ".txt"
        final_orbital_parameters_file = final_orbital_parameters_file_pw + str(no) + ".txt"
        
        f = open(composition_output_file, 'r')
        raw_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
        compositions = organize_compositions(raw_compositions) #organzies the data from the file            
        f.close()
        
        f = open(final_orbital_parameters_file, 'r')
        final_orbital_parameters = [line.split() for line in f.readlines()]
        fin_a = [float(final_orbital_parameters[i][2]) for i in range(len(final_orbital_parameters))]
        fin_e = [float(final_orbital_parameters[i][3]) for i in range(len(final_orbital_parameters))]
        f.close()
        
        for i, obj1 in enumerate(compositions):
            for j, obj2  in enumerate(final_orbital_parameters):
                if float(obj1[0]) == float(obj2[0]): # Makes sure the object is in both files
                    if no < 11:
                        distr_final_masses[0].append(float(obj1[1])*mass_conversion)
                        distr_final_cmfs[0].append(float(obj1[2]))
                        distr_final_sa[0].append(float(obj2[2]))
                    elif 11 <= no < 21:
                        distr_final_masses[1].append(float(obj1[1])*mass_conversion)
                        distr_final_cmfs[1].append(float(obj1[2]))
                        distr_final_sa[1].append(float(obj2[2]))
                    elif 21 <= no < 31:
                        distr_final_masses[2].append(float(obj1[1])*mass_conversion)
                        distr_final_cmfs[2].append(float(obj1[2]))
                        distr_final_sa[2].append(float(obj2[2]))
                    elif 31 <= no < 41:
                        distr_final_masses[3].append(float(obj1[1])*mass_conversion)
                        distr_final_cmfs[3].append(float(obj1[2]))
                        distr_final_sa[3].append(float(obj2[2]))
                    elif 41 <= no < 51:
                        distr_final_masses[4].append(float(obj1[1])*mass_conversion)
                        distr_final_cmfs[4].append(float(obj1[2]))
                        distr_final_sa[4].append(float(obj2[2]))
                    
    final_masses.append(distr_final_masses)
    final_cmfs.append(distr_final_cmfs)
    final_sa.append(distr_final_sa)
 
combined_data_points = [[] for i in range(len(distribution_list))]
for i in range(len(distribution_list)):
    for j in range(len(ef_markers)):
        for k in range(len(final_masses[i][j])):
            combined_data_point = [final_masses[i][j][k], final_cmfs[i][j][k], final_sa[i][j][k], ef_markers[j], ef_marker_sizes[j]] 
            combined_data_points[i].append(combined_data_point)
 

#for i in range(len(distribution_list)):
    #random.Random(1).shuffle(combined_data_points[i])
 
# Figure 20
cmf_min = 0.0
cmf_max = 1.0
point_opacity = 0.9
point_size = 15.0
#color_map = cmr.get_sub_cmap('cmr.sapphire_r', 0.2, 0.75)
color_map = 'viridis_r'
marker_type = 'o'
edge_color = 'none'
edge_width = 0.01

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, gridspec_kw={'width_ratios': [10, 10, 10.5]}, figsize=(13,5))
axes = [ax1, ax2, ax3]
fig.subplots_adjust(wspace=0.03)
ax1_legend_elements = [Line2D([], [], lw=0.0, color='black', marker=ef_markers[0], markersize=5.0, label='ef=3'),
                       Line2D([], [], lw=0.0, color='black', marker=ef_markers[1], markersize=5.0, label='ef=5'),
                       Line2D([], [], lw=0.0, color='black', marker=ef_markers[2], markersize=5.0, label='ef=7'),
                       Line2D([], [], lw=0.0, color='black', marker=ef_markers[3], markersize=5.0, label='ef=10'),
                       Line2D([], [], lw=0.0, color='black', marker=ef_markers[4], markersize=5.0, label='ef=15'),
                       Line2D([], [], color='black', linestyle=':', label='Embryo Mass', alpha=0.4),
                       Line2D([], [], color='black', linestyle='-.', label='Planetesimal Mass', alpha=0.4)] 
for i in range(len(combined_data_points[0])):
    p1 = ax1.scatter(combined_data_points[0][i][2], combined_data_points[0][i][0], c=combined_data_points[0][i][1], marker=combined_data_points[0][i][3], s=combined_data_points[0][i][4], edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)
    p2 = ax2.scatter(combined_data_points[1][i][2], combined_data_points[1][i][0], c=combined_data_points[1][i][1], marker=combined_data_points[1][i][3], s=combined_data_points[1][i][4], edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)
    p3 = ax3.scatter(combined_data_points[2][i][2], combined_data_points[2][i][0], c=combined_data_points[2][i][1], marker=combined_data_points[2][i][3], s=combined_data_points[2][i][4], edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)
    
    #p2 = ax2.scatter(final_sa[1][i], final_masses[1][i], c=final_cmfs[1][i], marker=marker, s=point_size, edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)
    #p3 = ax3.scatter(final_sa[2][i], final_masses[2][i], c=final_cmfs[2][i], marker=marker, s=point_size, edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)
    
"""for i, marker in enumerate(ef_markers):  
    p1 = ax1.scatter(final_sa[0][i], final_masses[0][i], c=final_cmfs[0][i], marker=marker, s=point_size, edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)
    p2 = ax2.scatter(final_sa[1][i], final_masses[1][i], c=final_cmfs[1][i], marker=marker, s=point_size, edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)
    p3 = ax3.scatter(final_sa[2][i], final_masses[2][i], c=final_cmfs[2][i], marker=marker, s=point_size, edgecolors=edge_color, linewidth=edge_width, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity, zorder=2)"""
for ax in axes:
    ax.axhline(.093, label = 'Embryo Mass', color = 'black', linestyle = ':', alpha=.7, zorder=1)
    ax.axhline(.0093, label = 'Planetesimal Mass', color = 'black', linestyle = '-.', alpha=.7, zorder=1)
    ax.set_xlabel("Semi-Major Axis (AU)", fontsize='x-large')
    ax.grid(alpha=0.7)
    ax.legend(handles=ax1_legend_elements, loc = 'upper right', framealpha=0.7, fontsize=6.5)
    
ax1.set_title('a)', fontsize='x-large')
ax2.set_title('b)', fontsize='x-large')
ax3.set_title('c)', fontsize='x-large')
ax1.set_ylabel('Mass ($M_{\u2295}$)',fontsize='x-large')
ax1.minorticks_on()
plt.xlim(0.0, 4.5)
plt.yscale('log')
cbar = fig.colorbar(p3, location='right', anchor=(0,0), pad=0.0, ticks=[0.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
cbar.set_label(label='CMF', size='x-large', labelpad=10.0)
cbar.minorticks_on()
plt.savefig("graphs/Final_Semi_Major_Axes.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/Final_Semi_Major_Axes.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/Final_Semi_Major_Axes.png', bbox_inches='tight', dpi=300)
    