"""
Produces Figure 18 from Ferich et al. (IN PREP)
Needs another version of the Organize Compositions 
function to extract the data
"""

import matplotlib.pyplot as plt
import cmasher as cmr
from astropy import units as u

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth
mass_expansion_factor = 1200 # Artificially inflates a particle's mass for graphing

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

composition_input_file = "DBCT_input/2step_DBCT_input.txt"

f = open(composition_input_file, 'r')
raw_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
init_compositions = organize_compositions(raw_compositions) #organzies the data from the file
init_masses = [obj[1]*mass_conversion*mass_expansion_factor for obj in init_compositions]
init_cmfs = [obj[2] for obj in init_compositions]
init_sa = [obj[3] for obj in init_compositions]
init_e = [obj[4] for obj in init_compositions]
f.close()

# Figure 18
plt.rcParams['axes.axisbelow'] = True # Makes sure grid is behind points
fig1, ax1 = plt.subplots(figsize=(6,5))
min_frac=0.1 
max_frac=0.7 
#color_map = cmr.sapphire_r
color_map = 'viridis_r'
plot1 = ax1.scatter(init_sa, init_e, marker='o', s=init_masses, c=init_cmfs, cmap=color_map, vmin=min_frac, vmax=max_frac)
ax1.set_xlabel('Semi-major Axis (AU)',fontsize='large')
ax1.set_ylabel("Eccentricity", fontsize='large')
ax1.set_xlim(0.0, 4.5)
ax1.set_ylim(-0.0001, 0.0101)
cbar1 = fig1.colorbar(plot1, location='right', anchor=(0,0.5), pad=0.0, ticks=[0.1, .2, .3, .4, .5, .6, .7])
cbar1.set_label(label='CMF', size='large')
cbar1.minorticks_on()
ax1.grid(alpha=0.7)
plt.savefig("graphs/2step_Initial_Disc.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/2step_Initial_Disc.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/2step_Initial_Disc.png', bbox_inches='tight', dpi=300)