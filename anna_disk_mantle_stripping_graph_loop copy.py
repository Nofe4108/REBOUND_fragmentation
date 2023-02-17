from astropy import units as u
import astropy.constants as constants
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import math
import time
import sys


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

final_hashes = []
final_masses = []
final_a = []
final_e = []
final_core_fracs = []
compositions = []

comp_output_file_pw = "mantle_stripping_output/low_exp_mantle_stripping_output"
final_orbital_parameters_file_pw = "final_orbital_parameters/final_orbital_parameters"
fig1_file = 'graphs/low_exp_all_final_core_fracs.pdf'
fig2_file = 'graphs/low_exp_all_final_planets.pdf'

file_range = np.arange(1,7,1)

for no in file_range:
    comp_output_file = comp_output_file_pw + str(no) + ".txt"
    final_orbital_parameters_file = final_orbital_parameters_file_pw + str(no) + ".txt"

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
                
    for c in comp:
        compositions.append(c)
    
    f = open(final_orbital_parameters_file, 'r')
    final_orbital_parameters = [line.split() for line in f.readlines()]
    final_hash = [int(final_orbital_parameters[i][0]) for i in range(len(final_orbital_parameters))]
    final_mass = [float(final_orbital_parameters[i][1])*334672.021419 for i in range(len(final_orbital_parameters))]
    fin_a = [float(final_orbital_parameters[i][2]) for i in range(len(final_orbital_parameters))]
    fin_e = [float(final_orbital_parameters[i][3]) for i in range(len(final_orbital_parameters))]
    for i in final_hash:
        final_hashes.append(i)
    for i in final_mass:
        final_masses.append(i)
    for i in fin_a:
        final_a.append(i)
    for i in fin_e:
        final_e.append(i)

        

final_core_fracs = []
final_masses_comp = []

for comp in compositions:
    final_masses_comp.append(comp[1]*334672.021419*1200)
    final_core_fracs.append(comp[2])



min_frac_final = min(final_core_fracs)
max_frac_final = max(final_core_fracs)

fig1, ax1 = plt.subplots()
ax1.scatter(final_masses, final_core_fracs, color = 'black')
ax1.axhline(0.3, label = 'Initial Core Fraction', color = 'tab:blue', linestyle = '--', alpha=.7)
ax1.axvline(.093, label = 'Initial Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax1.axvline(.0093, label = 'Initial Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax1.set_xlabel('Mass ($M_{\u2295}$)',fontsize='large')
ax1.set_ylabel("Core Fraction", fontsize='large')
plt.xscale('log')
ax1.set_ylim([0, 1.1])
plt.yticks(np.arange(0, 1.1, .1))
ax1.minorticks_on()
plt.grid()
plt.legend()

plt.savefig(fig1_file, bbox_inches='tight', pad_inches=0.25)

fig2, ax2 = plt.subplots(figsize=(8,5))
color_map = plt.get_cmap('jet_r')
plot2 = ax2.scatter(final_a, final_e, marker='o', s=final_masses_comp, c=final_core_fracs, cmap=color_map, vmin=min_frac_final, vmax=max_frac_final)
ax2.set_xlabel('Semi-major Axis (AU)',fontsize='large')
ax2.set_ylabel("Eccentricity", fontsize='large')
ax2.set_xlim(0.0, 4.5)
ax2.set_ylim(-0.01, 0.40)
cbar2 = fig2.colorbar(plot2, location='right', anchor=(0,0.5), pad=0.0)
cbar2.set_label(label='CMF', size='large')
cbar2.minorticks_on()
plt.grid()
ax2_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*1000)),
                          Line2D([], [], color='black', marker='o', lw=0.0, label='0.5 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(0.5*1000))]  
ax2_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*1200))]  
ax2_legend = plt.legend(handles=ax2_legend_elements, loc = 'upper right', prop={"size": 10})

plt.savefig(fig2_file, bbox_inches='tight', pad_inches=0.25)






