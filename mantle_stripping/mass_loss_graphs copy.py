#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 20:31:20 2023

@author: nferich
"""
from astropy import units as u
import astropy.constants as constants
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import matplotlib as matplotlib
import numpy as np
import math
import sys

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

#FILE PATHWAYS
ejection_file_pw = "ejections/ejections"
collision_report_file_pw = "new_collision_reports/new_collision_report"
comp_input_file_pw = "mantle_stripping_input/uni_mantle_stripping_input"
comp_output_file_pw = "mantle_stripping_output/uni_mantle_stripping_output"
ejection_compositions_pw = "ejections/uni_ejection_compositions" 
large_obj_collision_compositions_pw = "large_obj_collision_compositions/uni_large_obj_collision_compositions"
fig3_file = 'graphs/all_so_ejections.pdf'
fig4_file = 'graphs/all_so_collisions50.pdf'

mass_conversion = 334672.021419
artificial_expansion_factor = 500

#NUMBER OF FILES
file_range = np.arange(1,51,1)

#LISTS FOR EXTRACTED DATA
ejection_hashes = []
ejection_masses = []
ejection_times = []
ejection_core_fracs = []
large_obj_collision_times = []
large_obj_collision_flags = [] #indicates what object the projectile collided with - 1 for sun, 2 for jupiter, 3 for saturn
large_obj_collision_hashes = []
large_obj_collision_masses = []
large_obj_collision_core_fracs = []
initial_hashes = []
initial_masses = []
initial_core_fracs = []
initial_a = []
final_masses = []
final_core_fracs = []
ejection_so_masses = []
ejection_so_core_fracs = []
ejection_so_a = []
ejection_so_times = []
collision_so_masses = []
collision_so_core_fracs = []
collision_so_a = []
collision_so_times = []


################### START OF LOOP ##############################
for no in file_range:
    print(no)
    
    #FILES 
    ejection_file = ejection_file_pw + str(no) + ".txt"
    collision_report_file = collision_report_file_pw + str(no) + ".txt" 
    comp_input_file = comp_input_file_pw + str(1) + ".txt"
    comp_output_file = comp_output_file_pw + str(no) + ".txt"
    ejection_compositions_file = ejection_compositions_pw + str(no) + ".txt"
    large_obj_collision_compositions_file = large_obj_collision_compositions_pw + str(no) + ".txt"
    
    #EJECTION HASHES, MASSES, and TIMES
    ef = open(ejection_file, 'r')
    ejection_raw = [line.split() for line in ef.readlines()]
    ejection_time = [float(ejection_raw[i][1])/1.0e6 for i in range(len(ejection_raw))]
    ejection_hash = [int(ejection_raw[i][2]) for i in range(len(ejection_raw))]
    ejection_mass = [float(ejection_raw[i][3])*mass_conversion*artificial_expansion_factor for i in range(len(ejection_raw))]
    ejection_times.append(ejection_time)
    ejection_hashes.append(ejection_hash)
    ejection_masses.append(ejection_mass)
    
    #EJECTION CMFs
    ecf = open(ejection_compositions_file, 'r')
    ejection_comp_raw = [line.split() for line in ecf.readlines()]
    ejection_core_frac = [float(ejection_comp_raw[i][2]) for i in range(len(ejection_comp_raw))]
    ejection_core_fracs.append(ejection_core_frac)
    
    #INITIAL AND FINAL MASSES AND CMFs
    init_hash = []
    init_mass = []
    init_core_frac = []
    init_a = []
    fin_mass = []
    fin_core_frac = []
    
    cif = open(comp_input_file, 'r')
    init_compositions = [line.split() for line in cif.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    init_comp = organize_compositions(init_compositions) #organzies the data from the file
    for line in init_comp: #loops through each particle in the array
        if line[2] > 1.0 or line[2] < 0.0: #makes sure CMF isn't negative or greater than one
            print ('ERROR: CMF does not a realistic value')
            sys.exit(1)             
            try: 
                init_hashes = [int(x[0]) for x in init_compositions] 
            except:
                init_hashes = [x[0].value for x in init_compositions]
        init_hash.append(int(line[0]))
        init_mass.append(float(line[1])*mass_conversion*artificial_expansion_factor)
        init_core_frac.append(float(line[2]))
        init_a.append(float(line[3]))
        
    cof = open(comp_output_file, 'r')
    init_compositions = [line.split() for line in cof.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    fin_comp = organize_compositions(init_compositions) #organzies the data from the file
    for line in fin_comp: #loops through each particle in the array
        if line[2] > 1.0 or line[2] < 0.0: #makes sure CMF isn't negative or greater than one
            print ('ERROR: CMF does not a realistic value')
            sys.exit(1)             
            try: 
                init_hashes = [int(x[0]) for x in init_compositions] 
            except:
                init_hashes = [x[0].value for x in init_compositions]
        fin_mass.append(float(line[1])*mass_conversion*artificial_expansion_factor)
        fin_core_frac.append(float(line[2]))
    
    initial_hashes.append(init_hash)
    initial_masses.append(init_mass)
    initial_core_fracs.append(init_core_frac)
    initial_a.append(init_a)
    final_masses.append(fin_mass)
    final_core_fracs.append(fin_core_frac)
    
    #LARGE OBJECT COLLISIONS
    loccf = open(large_obj_collision_compositions_file, 'r')
    large_obj_collision_compositions_raw = [line.split() for line in loccf.readlines()]
    large_obj_collision_time = [float(large_obj_collision_compositions_raw[i][0])/1.0e6 for i in range(len(large_obj_collision_compositions_raw))]
    large_obj_collision_hash = [float(large_obj_collision_compositions_raw[i][2]) for i in range(len(large_obj_collision_compositions_raw))]
    large_obj_collision_mass = [float(large_obj_collision_compositions_raw[i][3])*mass_conversion*artificial_expansion_factor for i in range(len(large_obj_collision_compositions_raw))]
    large_obj_collision_core_frac = [float(large_obj_collision_compositions_raw[i][4]) for i in range(len(large_obj_collision_compositions_raw))]
    
    large_obj_collision_times.append(large_obj_collision_time)
    large_obj_collision_hashes.append(large_obj_collision_hash)
    large_obj_collision_masses.append(large_obj_collision_mass)
    large_obj_collision_core_fracs.append(large_obj_collision_core_frac)
    for collision in large_obj_collision_compositions_raw:
        if float(collision[1]) >= 0.98:
            large_obj_collision_flags.append(1)
        elif 1.0 > float(collision[1]) > 5.0e-4:
            large_obj_collision_flags.append(2)
        elif 5.0e-4 > float(collision[1]) > 1.0e-4:
            large_obj_collision_flags.append(3)
    
            
    #EJECTED STARTING MARTERIAL
    ejec_so_mass = []
    ejec_so_core_frac = []
    ejec_so_a = []
    ejec_so_time = []
    
    for i,hsh1 in enumerate(init_hash):
        for j,hsh2 in enumerate(ejection_hash):
            if int(hsh1) == int(hsh2): #finds starting object that was ejected
                if init_core_frac[i] == ejection_core_frac[j] and init_mass[i] == ejection_mass[j]: #if object is ejected without undergoing collision
                   ejec_so_mass.append(ejection_mass[j])  
                   ejec_so_core_frac.append(ejection_core_frac[j]) 
                   ejec_so_a.append(init_a[i]) 
                   ejec_so_time.append(ejection_time[j])
                   
    ejection_so_masses.append(ejec_so_mass)
    ejection_so_core_fracs.append(ejec_so_core_frac)
    ejection_so_a.append(ejec_so_a)
    ejection_so_times.append(ejec_so_time)
    
    #COLLISION STARTING MARTERIAL
    coll_so_mass = []
    coll_so_core_frac = []
    coll_so_a = []
    coll_so_time = []
    for i,hsh1 in enumerate(init_hash):
        for j,hsh2 in enumerate(large_obj_collision_hash):
            if hsh1 == hsh2: #finds starting object that was ejected
                if init_core_frac[i] == large_obj_collision_core_frac[j] and init_mass[i] == large_obj_collision_mass[j]: #if object is ejected without undergoing collision
                   coll_so_mass.append(large_obj_collision_mass[j])  
                   coll_so_core_frac.append(large_obj_collision_core_frac[j]) 
                   coll_so_a.append(init_a[i]) 
                   coll_so_time.append(large_obj_collision_time[j])
            
    collision_so_masses.append(coll_so_mass)
    collision_so_core_fracs.append(coll_so_core_frac)
    collision_so_a.append(coll_so_a)
    collision_so_times.append(coll_so_time)
    
    #CLOSE FILES
    ef.close()
    cif.close()
    cof.close()
    ecf.close()
    loccf.close()
############ END OF LOOP ########################


############## FINDING TOTALS ####################
initial_total_masses = [sum(initial_masses[i])/artificial_expansion_factor for i in range(len(initial_masses))]  
initial_total_core_masses = [sum(initial_masses[i])*sum(initial_core_fracs[i])/(len(initial_core_fracs[i])*artificial_expansion_factor) for i in range(len(initial_masses))] 
initial_avg_core_fracs = [sum(initial_core_fracs[i])/len(initial_core_fracs[i]) for i in range(len(initial_masses))] 
final_total_masses = [sum(final_masses[i])/artificial_expansion_factor for i in range(len(final_masses))]  
final_total_core_masses = [sum(final_masses[i])*sum(final_core_fracs[i])/(len(final_core_fracs[i])*artificial_expansion_factor) for i in range(len(final_masses))] 
total_ejected_masses = [sum(ejection_masses[i])/artificial_expansion_factor for i in range(len(ejection_masses))]   
total_mass_differences = [initial_total_masses[i]-final_total_masses[i] for i in range(len(initial_total_core_masses))]
total_collision_mass = [sum(large_obj_collision_masses[i])/artificial_expansion_factor for i in range(len(large_obj_collision_masses))]


################ GRAPHING #######################
min_frac_ejection = min(ejection_core_fracs[0])
max_frac_ejection = max(ejection_core_fracs[0])

################ FIGURE 1 #######################
"""fig1, ax1 = plt.subplots(figsize=(12,10))
color_map = plt.get_cmap('jet_r')
plot1 = ax1.scatter(ejection_times, ejection_core_fracs, marker='o')
ax1.set_xlabel('Collision time (yr)',fontsize='large')
ax1.set_ylabel("EMass(Me)", fontsize='large')
ax1.set_xlim(0.0, 0.2e7)
#ax1.set_ylim(-0.01, 0.40)
plt.grid()
ax1_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*1000)),
                          Line2D([], [], color='black', marker='o', lw=0.0, label='0.5 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(0.5*1000))]  
ax1_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*1200))]  
ax1_legend = plt.legend(handles=ax1_legend_elements, loc = 'upper right', prop={"size": 10})

########## FIGURE 2 ###############
fig2, ax2 = plt.subplots(figsize=(12,10))
plot2 = ax2.scatter(large_obj_collision_masses[0], large_obj_collision_core_fracs[0], marker='o')
ax2.set_xlabel('Collision time (yr)',fontsize='large')
ax2.set_ylabel("EMass(Me)", fontsize='large')
#ax1.set_xlim(0.0, 4.5)
#ax1.set_ylim(-0.01, 0.40)
plt.grid()"""

########## FLATTENED LISTS ###############
all_ejection_so_masses = [mass for sublist in ejection_so_masses for mass in sublist]
all_ejection_so_core_fracs = [core_frac for sublist in ejection_so_core_fracs for core_frac in sublist]
all_ejection_so_a = [a for sublist in ejection_so_a for a in sublist]
all_ejection_so_times = [time for sublist in ejection_so_times for time in sublist]
all_collision_so_masses = [mass for sublist in collision_so_masses for mass in sublist]
all_collision_so_core_fracs = [core_frac for sublist in collision_so_core_fracs for core_frac in sublist]
all_collision_so_a = [a for sublist in collision_so_a for a in sublist]
all_collision_so_times = [time for sublist in collision_so_times for time in sublist]

all_ejection_so_a_embryo = []
all_ejection_so_a_planetesimal = []

for i in range(len(all_collision_so_masses)):
    if all_collision_so_masses[i] < 0.09*artificial_expansion_factor:
        all_ejection_so_a_planetesimal.append(all_collision_so_a[i])
    else:
        all_ejection_so_a_embryo.append(all_collision_so_a[i])
        
all_ejection_so_a_nested  = [all_ejection_so_a_planetesimal, all_ejection_so_a_embryo]
print(all_ejection_so_masses[0:5])
########## FIGURE 3 ###############
bp_colors = ['tab:blue', 'tab:orange']
hist_type = 'barstacked'
n_bins = np.linspace(.5, 4.0, num=36, endpoint=True)
fig3, ax3 = plt.subplots(figsize=(6,5))
plot3 = ax3.hist(all_ejection_so_a_nested, n_bins, histtype=hist_type, color=bp_colors, lw=3.0, alpha = 1.0)
ax3.set_xlabel('Initial Semi-major Axis (AU)',fontsize='large')
ax3.set_ylabel("Number of Ejections", fontsize='large')
#ax3.set_xlim(-0.01, 2.01)
#ax3.set_ylim(-0.01, 0.40)
ax3.minorticks_on()
plt.grid()
ax3_legend_elements = [matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:blue', facecolor='tab:blue', label='planetesimal'),
                       matplotlib.patches.Rectangle((0,0), width=0.1, height=0.1, edgecolor='tab:orange', facecolor='tab:orange', label='embryo')]
legend = plt.legend(handles=ax3_legend_elements, loc = 'upper left', framealpha = .7)
plt.gca().add_artist(legend)
plt.savefig(fig3_file, bbox_inches='tight', pad_inches=0.01)


"""########## FIGURE 3 ###############
fig3, ax3 = plt.subplots(figsize=(7,5))
plot3 = ax3.scatter(all_ejection_so_times, all_ejection_so_a, s=all_ejection_so_masses, color='black', marker='o')
ax3.set_xlabel('Time of Ejection (Myr)',fontsize='large')
ax3.set_ylabel("Initial Semi-major Axis (AU)", fontsize='large')
ax3_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*artifical_expansion_factor))]                                                           
ax3_legend = plt.legend(handles=ax3_legend_elements, loc = 'upper right', prop={"size": 10})
ax3.set_xlim(-0.01, 2.01)
#ax3.set_ylim(-0.01, 0.40)
plt.grid()
plt.savefig(fig3_file, bbox_inches='tight', pad_inches=0.25, dpi=250)"""



"""########## FIGURE 4 ###############
fig4, ax4 = plt.subplots(figsize=(11,8))
plot4 = ax4.scatter(all_collision_so_times, all_collision_so_a, s=all_collision_so_masses, marker='o')
plt.title('Sun/Jupiter/Saturn Collisions', fontsize='large')
ax4.set_xlabel('Time of Collision (Myr)',fontsize='large')
ax4.set_ylabel("Initial Semi-major Axis (AU)", fontsize='large')
ax4_legend_elements = [Line2D([], [], color='black', marker='o', lw=0.0, label='0.1 $M_{\u2295}$', markerfacecolor='black', markersize=np.sqrt(.1*1200))]                                                           
ax4_legend = plt.legend(handles=ax3_legend_elements, loc = 'upper right', prop={"size": 10})
#ax3.set_xlim(0.3, 4.2)
#ax3.set_ylim(-0.01, 0.40)
plt.grid()
plt.savefig(fig4_file, bbox_inches='tight', pad_inches=0.25, dpi=250)"""
