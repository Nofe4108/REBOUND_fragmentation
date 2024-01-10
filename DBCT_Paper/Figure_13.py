"""
Produces Figure 13 from Ferich et al. (IN PREP)
Needs another version of the Organize Compositions 
function to extract the data
"""

import matplotlib.pyplot as plt
import cmasher as cmr
from astropy import units as u

def organize_compositions(init_compositions):
    """Data organizing function
    
    Takes in raw string data extracted from mantle stripping input file
    Creates a list filled with properties of each particle
    The properties include hash, mass, and CMF
    
    Parameters:
    init_compositions (list) -- raw particle data from input file

    Returns:
    compositions (list) -- nested list with properly formatted particle data
    """
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

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
mass_expansion_factor = 1200 # Artificially inflates a particle's mass for graphing
distribution_list = ['3step', 'lin', 'exp'] #Used to create full pathways to input files

# Lists used for grqphing
initial_masses = []
initial_cmfs = []
initial_sa = []
initial_e = []

# File pathways
composition_input_file1 = "DBCT_input/3step_DBCT_input.txt"
composition_input_file2 = "DBCT_input/lin_DBCT_input.txt"
composition_input_file3 = "DBCT_input/exp_DBCT_input.txt"
DBCT_input_file_pw = "DBCT_input/"
DBCT_input_file_name = "_DBCT_input"

for distr in distribution_list:
    composition_input_file = DBCT_input_file_pw + distr + DBCT_input_file_name + ".txt"
    f = open(composition_input_file, 'r')
    raw_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    init_compositions = organize_compositions(raw_compositions) #organzies the data from the file
    init_masses = [obj[1]*mass_conversion*mass_expansion_factor for obj in init_compositions]
    init_cmfs = [obj[2] for obj in init_compositions]
    init_sa = [obj[3] for obj in init_compositions]
    init_e = [obj[4] for obj in init_compositions]
    initial_masses.append(init_masses)
    initial_cmfs.append(init_cmfs)
    initial_sa.append(init_sa)
    initial_e.append(init_e)
    f.close()

min_initial_cmf = min(initial_cmfs[2]) # Found in the exponential distriubtion
max_initial_cmf = max(initial_cmfs[2]) # Found in the exponential distriubtion

# Figure 13
color_map = cmr.sapphire_r
fig1, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, gridspec_kw={'width_ratios': [10, 10, 10.5]}, figsize=(10,5))
fig1.subplots_adjust(wspace=0)
p1 = ax1.scatter(initial_sa[0], initial_e[0], marker='o', s=initial_masses[0], c=initial_cmfs[0], cmap=color_map, vmin=min_initial_cmf, vmax=max_initial_cmf+.01)
p2 = ax2.scatter(initial_sa[1], initial_e[1], marker='o', s=initial_masses[1], c=initial_cmfs[1], cmap=color_map, vmin=min_initial_cmf, vmax=max_initial_cmf+.01)
p3 = ax3.scatter(initial_sa[2], initial_e[2], marker='o', s=initial_masses[2], c=initial_cmfs[2], cmap=color_map, vmin=min_initial_cmf, vmax=max_initial_cmf+.01)
ax1.minorticks_on()
ax1.set_title('a)', fontsize='x-large')
ax2.set_title('b)', fontsize='x-large')
ax3.set_title('c)', fontsize='x-large')
fig1.supxlabel('Semi-major Axis (AU)', fontsize='x-large')
ax1.set_ylabel("Eccentricity", fontsize='x-large', labelpad=10.0)
plt.xlim(0.0, 4.5)
plt.ylim(0.000, 0.0101)
cbar = fig1.colorbar(p3, location='right', anchor=(0,0), pad=0.0, ticks=[.1, .2, .3, .4, .5, .6, .7])
cbar.set_label(label='CMF', size='x-large', labelpad=10.0)
cbar.minorticks_on()
plt.savefig("graphs/fig13.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig13.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig13.png', bbox_inches='tight', dpi=300)
