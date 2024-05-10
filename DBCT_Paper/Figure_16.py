"""
Produces Figure 16 from Ferich et al. (IN PREP)
Needs to use a different form of the Organize
Compositions function from the DBCT to 
extract the necessary data
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy import units as u
import numpy as np

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

mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth

# ND stands for non-disruptive - this list will include objects that did not undergo a mass or composition changing collision
nd_ejection_masses = [[], []] # First inner list is for planetesimals (Tab:Blue) and second list is for embryos (tab:orange)
nd_ejection_sa = [[], []] # Original emi-major axes of objects when they get ejeted - first list is for planetesimals, second if for embryos

# File Pathways
collision_report_file_pw = "new_collision_reports/new_collision_report"
composition_input_file_pw = "DBCT_input/uni_DBCT_input"
ejection_file_pw = "ejections/ejections"

file_range = np.arange(1,51,1) # Number of simulations the code produces data for

for no in file_range:    
    
    nd_ejec_masses = []
    nd_ejec_sa = []
    
    ejection_file = ejection_file_pw + str(no) + ".txt"
    collision_report_file = collision_report_file_pw + str(no) + ".txt" 
    composition_input_file = composition_input_file_pw + ".txt"
    
    # Extracts compositional and orbital data of objects in the initial disc
    f = open(composition_input_file, 'r')
    raw_compositions = [line.split() for line in f.readlines()]
    compositions = organize_compositions(raw_compositions) # This list keeps track of all the particle data
    init_hashes = [obj[0] for obj in compositions]
    init_masses = [obj[1]*mass_conversion for obj in compositions]
    init_sa = [obj[3] for obj in compositions]
    f.close()
    
    # Extracts the masses and hashes of ejected objects
    f = open(ejection_file, 'r')
    ejection_raw = [line.split() for line in f.readlines()]
    ejection_hashes = [int(ejec[2]) for ejec in ejection_raw]
    #print(len(ejection_hashes))
    ejection_masses = [float(ejec[3])*mass_conversion for ejec in ejection_raw]
    f.close()
    
    # Extracts collisional data
    f = open(collision_report_file, 'r')
    collision_blocks = f.read().split("\n")
    collisions = [block.split() for block in collision_blocks if len(block) > 0]
    collision_types = [int(collision[1]) for collision in collisions]
    targ_hashes = [int(collision[2]) for collision in collisions]
    proj_hashes = [int(collision[4]) for collision in collisions]
    f.close()
    
    # Loop used to determine
    for i,ejec_hsh in enumerate(ejection_hashes):
        dc_flag = 0 # Flag will indicate if an object undergoes a disruptive or mass changing collision
        for j,init_hsh in enumerate(init_hashes):
            if ejec_hsh == init_hsh: # Checks to see if the object was in the initial disc
                for k, collision_type in enumerate(collision_types):
                    if collision_type == 1 or collision_type == 2 or collision_type == 3 or collision_type == 4:
                        if ejec_hsh == targ_hashes[k] or ejec_hsh == proj_hashes[k]: # Checks if an object was the target or projectile of composition-changing collision
                            dc_flag += 1
                            break
                if dc_flag == 0:
                    nd_ejec_masses.append(ejection_masses[i])
                    nd_ejec_sa.append(init_sa[j])     
   
    for i, mass in enumerate(nd_ejec_masses):
        if mass < 0.01: # If object is a planetesimal
            nd_ejection_masses[0].append("tab:blue")
            nd_ejection_sa[0].append(nd_ejec_sa[i])
        else: # If object is an embryo
            nd_ejection_masses[1].append("tab:orange")
            nd_ejection_sa[1].append(nd_ejec_sa[i])
                 
# Figure 16
bp_colors = ['tab:blue', 'tab:orange']
hist_type = 'barstacked'
n_bins = np.linspace(.5, 4.0, num=36, endpoint=True)
fig1, ax1 = plt.subplots(figsize=(6,5))
plot1 = ax1.hist(nd_ejection_sa, n_bins, histtype=hist_type, color=bp_colors, lw=0.1, alpha = 1.0, edgecolor='black')
ax1.set_xlabel('Initial Semi-major Axis (AU)',fontsize='large')
ax1.set_ylabel("Number of Ejections", fontsize='large')
ax1.minorticks_on()
ax1.grid(alpha=0.7)
ax1_legend_elements = [mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:blue', facecolor='tab:blue', label='planetesimal'),
                       mpatches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:orange', facecolor='tab:orange', label='embryo')]
legend = plt.legend(handles=ax1_legend_elements, loc = 'upper left', framealpha = .7)
plt.savefig("graphs/fig16.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig16.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig16.png', bbox_inches='tight', dpi=300)
           
            