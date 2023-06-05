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

comp_output_file_pw = "comp_tracking_output/"
comp_output_file_name = "_comp_tracking_output"
fig_file = 'graphs/comp_tracking_all_final_core_fracs.pdf'


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
        
                

############# END OF LOOP ############################        

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, gridspec_kw={'width_ratios': [10, 10, 10]}, figsize=(13,5))
fig.subplots_adjust(wspace=0.03)
ax1.minorticks_on()

ax1.scatter(final_masses[0], final_core_fracs[0], s=3.0, color = 'black')
ax2.scatter(final_masses[1], final_core_fracs[1], s=3.0, color = 'black')
ax3.scatter(final_masses[2], final_core_fracs[2], s=3.0, color = 'black')
ax1.axvline(.093, label = 'Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax1.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax2.axvline(.093, label = 'Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax2.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax3.axvline(.093, label = 'Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax3.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax1.set_title('a)', fontsize='x-large')
ax2.set_title('b)', fontsize='x-large')
ax3.set_title('c)', fontsize='x-large')
ax1.set_xlabel('Mass ($M_{\u2295}$)', fontsize='x-large', labelpad=10.0)
ax2.set_xlabel('Mass ($M_{\u2295}$)', fontsize='x-large', labelpad=10.0)
ax3.set_xlabel('Mass ($M_{\u2295}$)', fontsize='x-large', labelpad=10.0)
ax1.set_ylabel("CMF", fontsize='x-large', labelpad=10.0)
plt.xscale('log')
plt.ylim(-0.01, 1.01)
plt.yticks(np.arange(0.0, 1.1, 0.1))
ax1.minorticks_on()
ax1.grid()
ax2.grid()
ax3.grid()
ax1.legend(framealpha=1.0, fontsize=8.0)
ax2.legend(framealpha=1.0, fontsize=8.0)
ax3.legend(framealpha=1.0, fontsize=8.0)

plt.savefig(fig_file, bbox_inches='tight', pad_inches=0.01)
