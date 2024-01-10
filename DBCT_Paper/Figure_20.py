"""
Produces Figure 20 from Ferich et al. (IN PREP)
"""

import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
from astropy import units as u
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
distribution_list = ['3step', 'lin', 'exp'] #Used to create full pathways to output files

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
    distr_final_masses = []
    distr_final_cmfs = []
    distr_final_sa = []
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
        
        if len(compositions) == len(final_orbital_parameters): # Ensures the orbital parameter and comp output file have same number of objects
            for i, obj in enumerate(compositions):
                distr_final_masses.append(float(obj[1])*mass_conversion)
                distr_final_cmfs.append(float(obj[2]))
                distr_final_sa.append(float(final_orbital_parameters[i][2]))
                
    final_masses.append(distr_final_masses)
    final_cmfs.append(distr_final_cmfs)
    final_sa.append(distr_final_sa)
                
 
# Figure 20
cmf_min = 0.0
cmf_max = 1.0
point_opacity = 0.9
point_size = 15.0
color_map = cmr.sapphire_r
marker_type = 'o'
edge_color = 'none'

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, gridspec_kw={'width_ratios': [10, 10, 10.5]}, figsize=(13,5))
axes = [ax1, ax2, ax3]
fig.subplots_adjust(wspace=0.03)
p1 = ax1.scatter(final_sa[0], final_masses[0], c=final_cmfs[0], marker=marker_type, s=point_size, edgecolors=edge_color, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity)
p2 = ax2.scatter(final_sa[1], final_masses[1], c=final_cmfs[1], marker=marker_type, s=point_size, edgecolors=edge_color, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity)
p3 = ax3.scatter(final_sa[2], final_masses[2], c=final_cmfs[2], marker=marker_type, s=point_size, edgecolors=edge_color, cmap=color_map, vmin=cmf_min, vmax=cmf_max, alpha=point_opacity)
for ax in axes:
    ax.axhline(.093, label = 'Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
    ax.axhline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
    ax.set_xlabel("Semi-Major Axis (AU)", fontsize='x-large')
    ax.grid(alpha=0.7)
    ax.legend(framealpha=1.0, fontsize=7.5)

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
plt.savefig("graphs/fig20.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig20.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig20.png', bbox_inches='tight', dpi=300)
    