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
from scipy.stats import norm
import statistics


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


def determine_ef(file_no, collision_index):
    if file_no <= 10:
        ef3_collision_types[collision_index] += 1
    elif 11 <= file_no <= 20:
        ef5_collision_types[collision_index] += 1
    elif 21 <= file_no <= 30:
        ef7_collision_types[collision_index] += 1
    elif 31 <= file_no <= 40:
        ef10_collision_types[collision_index] += 1
    elif 41 <= file_no <= 50:
        ef15_collision_types[collision_index] += 1
    

###### MAIN FUNCTION #########

final_hashes = []
final_masses = []
final_a = []
final_e = []
final_core_fracs = []
compositions = []
expansion_factors = []
ef_colors = []
bp_final_masses = []
bp_CMF_exp3 = []
bp_CMF_exp5 = []
bp_CMF_exp7 = []
bp_CMF_exp10 = []
bp_CMF_exp15 = []
ef3_collision_types = np.array([0,0,0,0,0])
ef5_collision_types = np.array([0,0,0,0,0])
ef7_collision_types = np.array([0,0,0,0,0])
ef10_collision_types = np.array([0,0,0,0,0])
ef15_collision_types = np.array([0,0,0,0,0])
collision_types = ['elastic bounce','merger','partial accretion','partial erosion','super catastrophic']

collision_report_file_pw = "new_collision_reports/new_collision_report"
comp_output_file_pw = "mantle_stripping_output/2step_mantle_stripping_output"
final_orbital_parameters_file_pw = "final_orbital_parameters/final_orbital_parameters"
fig1_file = 'graphs/2step_all_final_core_fracs.pdf'
fig2_file = 'graphs/2step_all_final_planets.pdf'
fig3_file = 'graphs/2step_all_final_core_fracs_efs.pdf'
fig4_file = 'graphs/2step_planet_CMF_histogram.pdf'
fig5_file = 'graphs/2step_planet_CMF_mass_histogram.pdf'
fig6_file = 'graphs/2step_planet_CMF_gaussians.pdf'
fig7_file = 'graphs/ef_collision_types_bar.pdf'

file_range = np.arange(1,51,1)

mass_conversion = 334672.021419
artificial_mass_expansion = 1
bp_mass = .093

def determine_ef(file_no, collision_index):
    if file_no <= 10:
        ef3_collision_types[collision_index] += 1
    elif 11 <= file_no <= 20:
        ef5_collision_types[collision_index] += 1
    elif 21 <= file_no <= 30:
        ef7_collision_types[collision_index] += 1
    elif 31 <= file_no <= 40:
        ef10_collision_types[collision_index] += 1
    elif 41 <= file_no <= 50:
        ef15_collision_types[collision_index] += 1

############## START OF LOOP #########################
for no in file_range:
    
    collision_report_file = collision_report_file_pw + str(no) + ".txt"    
    comp_output_file = comp_output_file_pw + str(no) + ".txt"
    final_orbital_parameters_file = final_orbital_parameters_file_pw + str(1) + ".txt"
    file = open(collision_report_file, 'r')
    blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
    blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list
    for i in range(len(blocks)): #iterates through each value in blocks list - THIS IS A VERY BIG LOOP THAT CONTAINS ALL OF THE MASS TRANSFER DECISION MAKING
        block = blocks[i].split() #seperates each long string full of the collision data into its own list to be parsed through
        time = float(block[0]) #time of collision is first value
        collision_type = int(block[1]) #type of collision is second
        if collision_type == 0:
            determine_ef(no, 0)
        elif collision_type == 1:
            determine_ef(no, 1)
        elif collision_type == 2:
            determine_ef(no, 2)
        elif collision_type == 3:
            determine_ef(no, 3)
        else:
            determine_ef(no, 4)
    
    
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
        if no < 11:
            expansion_factors.append('ef=3')
            ef_colors.append('tab:blue')
            if c[1]*mass_conversion*artificial_mass_expansion > bp_mass:
                bp_CMF_exp3.append(c[2])
        elif 11 <= no < 21:
            expansion_factors.append('ef=5')
            ef_colors.append('tab:orange')
            if c[1]*mass_conversion*artificial_mass_expansion > bp_mass:
                bp_CMF_exp5.append(c[2])
        elif 21 <= no < 31:
            expansion_factors.append('ef=5')
            ef_colors.append('tab:green')
            if c[1]*mass_conversion*artificial_mass_expansion > bp_mass:
                bp_CMF_exp7.append(c[2])
        elif 31 <= no < 41:
            expansion_factors.append('ef=10')
            ef_colors.append('tab:purple')
            if c[1]*mass_conversion*artificial_mass_expansion > bp_mass:
                bp_CMF_exp10.append(c[2])
        elif 41 <= no < 51:
            expansion_factors.append('ef=15')
            ef_colors.append('tab:red')
            if c[1]*mass_conversion*artificial_mass_expansion > bp_mass:
                bp_CMF_exp15.append(c[2])
    
    f = open(final_orbital_parameters_file, 'r')
    final_orbital_parameters = [line.split() for line in f.readlines()]
    final_hash = [int(final_orbital_parameters[i][0]) for i in range(len(final_orbital_parameters))]
    final_mass = [float(final_orbital_parameters[i][1])*mass_conversion for i in range(len(final_orbital_parameters))]
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

############# END OF LOOP ############################        

final_core_fracs = []
final_masses_comp = []

for comp in compositions:
    final_masses_comp.append(comp[1]*mass_conversion*artificial_mass_expansion)
    final_core_fracs.append(comp[2])
    if comp[1]*mass_conversion*artificial_mass_expansion > bp_mass:
        bp_final_masses.append(comp[1]*mass_conversion*artificial_mass_expansion)

bp_CMFs = [bp_CMF_exp3, bp_CMF_exp5, bp_CMF_exp7, bp_CMF_exp10, bp_CMF_exp15]
bp_CMFs_flat = [i for lst in bp_CMFs for i in lst]

min_frac_final = min(final_core_fracs)
max_frac_final = max(final_core_fracs)

fig1, ax1 = plt.subplots(figsize=(6,5))
ax1.scatter(final_masses_comp, final_core_fracs, s=5.0, color = 'black')
#ax1.axhline(0.3, label = 'Initial CMF', color = 'tab:blue', linestyle = '--', alpha=.7)
ax1.axvline(.093, label = 'Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax1.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax1.set_xlabel('Mass ($M_{\u2295}$)',fontsize='large')
ax1.set_ylabel("CMF", fontsize='large')
plt.xscale('log')
ax1.set_ylim([0, 1.1])
plt.yticks(np.arange(0, 1.1, .1))
ax1.minorticks_on()
plt.grid()
plt.legend()

plt.savefig(fig1_file, bbox_inches='tight', pad_inches=0.01)
"""
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


fig3, ax3 = plt.subplots(figsize=(13,10))


ax3.scatter(final_masses_comp, final_core_fracs, s=20.0, color = ef_colors)
ax3.axhline(0.3, label = 'Initial Core Fraction', color = 'tab:blue', linestyle = '--', alpha=.7)
ax3.axvline(.093, label = 'Initial Embryo Mass', color = 'tab:orange', linestyle = '--', alpha=.7)
ax3.axvline(.0093, label = 'Initial Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.7)
ax3.set_xlabel('Mass ($M_{\u2295}$)',fontsize='large')
ax3.set_ylabel("Core Fraction", fontsize='large')
plt.xscale('log')
ax3.set_ylim([0, 1.1])
plt.yticks(np.arange(0, 1.1, .1))
ax3.minorticks_on()
plt.grid()

ax3_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='ef=3', markerfacecolor='tab:blue', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='ef=5', markerfacecolor='tab:orange', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:green', marker='o', lw=0.0, label='ef=7', markerfacecolor='tab:green', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:purple', marker='o', lw=0.0, label='ef=10', markerfacecolor='tab:purple', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:red', marker='o', lw=0.0, label='ef=15', markerfacecolor='tab:red', markersize=np.sqrt(.1*1000)),]
second_legend = plt.legend(handles=ax3_legend_elements, loc = 'upper right', framealpha = .7)
plt.gca().add_artist(second_legend)
ax3.legend(loc=4, framealpha=0.6) 

plt.savefig(fig3_file, bbox_inches='tight', pad_inches=0.25)


fig4, ax4 = plt.subplots(figsize=(6,5))

n_bins = np.linspace(.25, .36, num=22, endpoint=False)
#n_bins = np.linspace(.05, .55, num=50, endpoint=False)
#n_bins = np.linspace(0.00, 0.55, num=55, endpoint=False)
hist_type = 'barstacked'
bp_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:red']



ax4.hist(bp_CMFs, bins=n_bins, histtype=hist_type, lw=3.0, color=bp_colors, alpha = 1.0)
ax4.set_xlim([0.24, 0.36])
ax4.set_xlabel('CMF',fontsize='large')
ax4.set_ylabel('Number of Planets (M>0.093 $M_{\u2295}$)', fontsize='large')
#plt.yticks(np.arange(0, 1.1, .1))
ax4.minorticks_on()
plt.grid()

ax4_legend_elements = [matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:blue', facecolor='tab:blue', label='ef=3'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:orange', facecolor='tab:orange', label='ef=5'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:green', facecolor='tab:green', label='ef=7'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:purple', facecolor='tab:purple', label='ef=10'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:red', facecolor='tab:red', label='ef=15')]
legend = plt.legend(handles=ax4_legend_elements, loc = 'upper right', framealpha = .7)
plt.gca().add_artist(legend)

plt.savefig(fig4_file, bbox_inches='tight', pad_inches=0.01)

fig5, ax5 = plt.subplots(figsize=(10,8))

n_bins = [20, 20]
#n_bins = [25, 20]
hist_type = 'barstacked'
bin_limits = [[.25, .35],[0.0, 2.0]]
#bin_limits = [[.05, .55],[0.0, 2.0]]

fig5 = ax5.hist2d(bp_CMFs_flat, bp_final_masses, bins=n_bins, range=bin_limits, cmap=plt.cm.jet)
ax5.set_xlim([0.25, 0.35])
ax5.set_ylim([0.1, 2.0])
ax5.set_xlabel('CMF (M>0.1 $M_{\u2295}$)',fontsize='large')
ax5.set_ylabel('$M_{\u2295}$', fontsize='large')
#plt.yticks(np.arange(0, 1.1, .1))
ax5.minorticks_on()
#plt.grid()

cbar5 = plt.colorbar(fig5[3], location='right', anchor=(0,0.5), pad=0.0)


plt.savefig(fig5_file, bbox_inches='tight', pad_inches=0.25)

for i in range(len(bp_CMFs[0])):
    if bp_CMFs[0][i] < 0.1:
        bp_CMFs[0].pop(i)
        break
    
cmf_ef_means = [np.mean(bp_CMFs[i]) for i in range(len(bp_CMFs))]
cmf_ef_sds = [np.std(bp_CMFs[i]) for i in range(len(bp_CMFs))]
print(np.mean(bp_CMFs[0]))
print(np.std(bp_CMFs[3]))

# Plot between -10 and 10 with .001 steps.
x_axis = np.arange(0.25, 0.35, 0.001)
efs = ['3','5','7','10','15']
fig6, ax6 = plt.subplots(figsize=(9,7))

for i in range(len(bp_CMFs)):
    ax6.plot(x_axis, norm.pdf(x_axis, cmf_ef_means[i], cmf_ef_sds[i]), lw=3.0, label = 'ef='+efs[i], c=bp_colors[i])
plt.grid()
plt.legend()
ax6.set_xlabel('CMF (M>0.1 $M_{\u2295}$)',fontsize='large')
ax6.set_ylabel('$M_{\u2295}$', fontsize='large')

plt.savefig(fig6_file, bbox_inches='tight', pad_inches=0.01)"""
print(len(ef5_collision_types))
fig7, ax7 = plt.subplots(figsize=(9,7))
ax7.minorticks_on()
plt.bar(collision_types, ef3_collision_types, color='tab:blue')
plt.bar(collision_types, ef5_collision_types, bottom=ef3_collision_types, color='tab:orange')
plt.bar(collision_types, ef7_collision_types, bottom=ef3_collision_types+ef5_collision_types, color='tab:green')
plt.bar(collision_types, ef10_collision_types, bottom=ef3_collision_types+ef5_collision_types+ef7_collision_types, color='tab:purple')
plt.bar(collision_types, ef15_collision_types, bottom=ef3_collision_types+ef5_collision_types+ef7_collision_types+ef10_collision_types, color='tab:red')


ax7_legend_elements = [matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:blue', facecolor='tab:blue', label='ef=3'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:orange', facecolor='tab:orange', label='ef=5'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:green', facecolor='tab:green', label='ef=7'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:purple', facecolor='tab:purple', label='ef=10'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:red', facecolor='tab:red', label='ef=15')]
legend = plt.legend(handles=ax7_legend_elements, loc = 'upper right', framealpha = .7)
plt.gca().add_artist(legend)

plt.savefig(fig7_file, bbox_inches='tight', pad_inches=0.01)
