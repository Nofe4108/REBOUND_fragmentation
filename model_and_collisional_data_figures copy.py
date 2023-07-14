#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 16:16:47 2023

@author: nferich
"""

import numpy as np
import math
import time
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#hello = [x.strip(' ') for x in hello]

#SIM CONVERSIONS
distance_conversion = 1.49597870691e11 #m to au
mass_conversion = 1.9885e30 #kg to Msun
time_conversion = 31557600 #sec to yrs
obj_mass_conversion = 334672.021419
G = 39.478

core_density = 7874.0*(distance_conversion**3/mass_conversion) #kg m^-3 - density of iron
mantle_density = 3000.0*(distance_conversion**3/mass_conversion) #kg m^3 -  value used in the simulation

efs = []
B_list = []
ideal_ejecta_CMF = []
actual_ejecta_CMF = []
collision_types = []
velocities = []
disruptive_collision_velocities = []
disruptive_collision_target_masses = []
disruptive_collision_proj_masses = []
disruptive_collision_lr_masses = []
disruptive_collision_lr_cmfs = []
disruptive_collision_target_cmfs = []
disruptive_collision_fragment_masses = []
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
            elif i < 4: #The other values will be floats - doesn't add initial semi-major axis or eccentricity to list
                particle_data.append(float(line[i]))
        compositions.append(particle_data) #add all the data about the particle to the compositions list
    return(compositions)

def calc_core_radius(mass, core_frac, core_density):
    core_mass = mass*core_frac
    core_radius = ((3*core_mass)/(4*np.pi*core_density))**(1/3)
    return(core_radius)

def calc_radius(mass, core_frac, mantle_density, core_radius):
    mantle_frac = 1-core_frac
    mantle_mass = mass*mantle_frac
    radius = (core_radius**3+((3*mantle_mass)/(4*np.pi*mantle_density)))**(1/3)
    return(radius)
      
#ratio = (2/np.pi)*math.atan((2-impact_param_over_Rtarg)*Q_over_Qstar/damping_factor)   
start_time = time.time() #Timer to see how long running this code takes


###### MAIN FUNCTION #########
def track_composition(collision_report_file, composition_input_file, impact_parameter_file, ejection_file, velocity_file, min_core_collision_frac, max_core_collision_frac, file_range_idx): #Main function that gives the compositions that will be outputted to the file
    f = open(composition_input_file, 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) #organzies the data from the file
    for line in compositions: #loops through each particle in the array
        if line[2] > 1.0 or line[2] < 0.0: #makes sure CMF isn't negative or greater than one
            print ('ERROR: CMF does not a realistic value')
            sys.exit(1)        
    try: 
        init_hashes = [int(x[0]) for x in init_compositions] 
    except:
        init_hashes = [x[0].value for x in init_compositions]
    
     
    file = open(collision_report_file, 'r')
    blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
    blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list
    
    #Gets the impact parameters for each collision
    b_file = open(impact_parameter_file, 'r')
    b = b_file.read() #reads each line in file and splits its value into an array - these impact parameters are divided by R_targ
    impact_parameters_raw = b.split("\n")
    impact_parameters = [x.strip('b/Rtarg:     ') for x in impact_parameters_raw]
    
    #GETS COLLISIONAL VELOCITIES FOR GRAPHING
    v_file = open(velocity_file, 'r')
    v = v_file.read() #reads each line in file and splits its value into an array - these impact parameters are divided by R_targ
    velocities_raw = v.split("\n")
    collision_velocities = [x.strip('Vimp/Vesc:     ') for x in velocities_raw]

    destroyed_object_hashes = [] #list of objects that get destroyed in a collision
    
    #START OF BIG LOOP
    for i in range(len(blocks)): #iterates through each value in blocks list - THIS IS A VERY BIG LOOP THAT CONTAINS ALL OF THE MASS TRANSFER DECISION MAKING
        block = blocks[i].split() #seperates each long string full of the collision data into its own list to be parsed through
        time = float(block[0]) #time of collision is first value
        collision_type = int(block[1]) #type of collision is second
        if collision_type == 0: #This is just an elastic bounce so it just continues to the next set of collision data in the blocks list
            velocities.append(float(collision_velocities[i]))
            continue
        target_hash = int(block[2]) #hash of the target
        largest_remnant_mass = float(block[3]) #mass of the target after collision
        proj_hash = int(block[4]) #hash of projectile

        
        #This section is here in case a projectile collides with the star or Jupiter - make this better
        big_object_collision = 1
        for i in range(len(compositions)):
            if target_hash==compositions[i][0]:
                big_object_collision+= -1
                break
        if big_object_collision==1:
            destroyed_object_hashes.append(proj_hash)
            continue   

        targ_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==target_hash][0] #this searches through the compositions array to see where the target hash matches up with a hash in the array - variable saves this index so it can go right to the correct row for the target in the compositions array during future oprations
        proj_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==proj_hash][0]  #same thing as above except for projectile
        target_mass = float(compositions[targ_idx][1]) #gets the mass of target before collision - uses the targ_index variable to go to the right row
        proj_mass = float(compositions[proj_idx][1]) #same thing as above except for projectile
        target_core_frac = compositions[targ_idx][2] #list of composition fractions for target before the collision - [-no.species:] -no.species means it goes back from the end of the list to the versy first composition fraction and the colon means it adds all values from the list past that into a list 
        proj_core_frac = compositions[proj_idx][2] #same thing as above except for projectile
        no_frags = int((len(block)-5)/2) #gives number of frags formed in collision - gets rid of 5 values from length which are basically about the collision type and target and projectile then divides in half since the list contains both the fragment hash and its mass
        frag_hashes = [int(block[i*2+3]) for i in range(1,no_frags+1)] #list of the hashes of the fragments - jumps from hash to hash using (i*2+3) - range length uses no_frags variable to get the correct amount of hashes
        frag_masses = [float(block[i*2+4]) for i in range(1,no_frags+1)] #same thing as above except for a list of the fragment masses 
        target_radius = calc_core_radius(target_mass, 1.0, mantle_density) #radius of target in sim before collision
        proj_radius = calc_core_radius(proj_mass, 1.0, mantle_density) #radius of projectile in sim before collision
        impact_parameter = float(impact_parameters[i])*target_radius
        
        
        #This sequence creates an estimate for the core radius of the target and projectile since the sim doesn't track that
        diff_target_core_radius = calc_core_radius(target_mass, target_core_frac, core_density)
        diff_proj_core_radius = calc_core_radius(proj_mass, proj_core_frac, core_density)
        diff_target_radius = calc_radius(target_mass, target_core_frac, mantle_density, diff_target_core_radius)
        diff_proj_radius = calc_radius(proj_mass, proj_core_frac, mantle_density, diff_proj_core_radius)
        target_radius_ratio = diff_target_core_radius/diff_target_radius
        proj_radius_ratio =  diff_proj_core_radius/diff_proj_radius
        target_core_radius = target_radius_ratio*target_radius
        proj_core_radius = proj_radius_ratio*proj_radius
        
                                               
 ######################## PERFECT MERGER ##########################
        if collision_type == 1: #perfect merger
            velocities.append(collision_velocities[i])
            compositions[targ_idx][2] = ((target_core_frac*target_mass)+(proj_core_frac*proj_mass))/largest_remnant_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile composition and mass
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
            
 
 ####################### PARTIAL ACCRETION ######################
   
        if collision_type == 2: #partial accretion
        
            mass_accreted = largest_remnant_mass-target_mass #change in mass of the target after the collision - this time mass is added to target
            ejecta_core_frac = 0 #fraction of the accreted mass that is composed of core material (will depend on impact parameter)

            if proj_core_radius == 0: #makes sure the core radius of the projectile isn't 0
                ejecta_core_frac = 0
            else:
                slope = (max_core_collision_frac-min_core_collision_frac)/(-2*proj_core_radius)
            
                max_impact_parameter = target_radius+proj_core_radius #if impact parameter is less than this than it will strip off more core
                min_impact_parameter = target_radius-proj_core_radius #if impact parameter is greater than this than it will strip off less core
            
                #DETERMINES THE FRACTION OF CORE AND MANTLE THAT GETS LOST FROM PROJECTILE
                if impact_parameter >= max_impact_parameter: #if the impact parameter is large (usually means oblique angle impact)
                    ejecta_core_frac = min_core_collision_frac #least amount of core is lost because cross-section of target won't intersect with core of projectile
                elif max_impact_parameter > impact_parameter > min_impact_parameter: #if it's between these extreme values
                    ejecta_core_frac = min_core_collision_frac+(slope*(impact_parameter-max_impact_parameter)) #it will strip off a certain amount of core according to this line equation
                elif impact_parameter <= min_impact_parameter: #if the impact parameter is small (usually means a more head-on impact)
                    ejecta_core_frac = max_core_collision_frac #most amount of core is lost because cross-section of target will fully intersect with core of projectile
            
            
            ideal_mass_accreted= [(mass_accreted*ejecta_core_frac),(mass_accreted*(1-ejecta_core_frac))] #first is core mass, second if mantle mass - how much of each layer would ideally be accreted if the projectile has the right composition
            proj_layer_mass = [proj_mass*proj_core_frac, proj_mass*(1-proj_core_frac)] #how much mass is in each layer of projectile - core first, then mantle
            
            mass_accreted_fractions = [] #will hold the actual fractions of each layer that get accreted from the projectile

            #these two statements make sure the correct amount of transferred mass is obtained from each layer - if there's not enough mass in a layer it takes the rest from the other layer
            if proj_layer_mass[0] >= ideal_mass_accreted[0]:
                mass_accreted_fractions.append(ideal_mass_accreted[0])
            else:
                mass_accreted_fractions.append(proj_layer_mass[0])
                ideal_mass_accreted[1] += (ideal_mass_accreted[0]-proj_layer_mass[0])
                
            if proj_layer_mass[1] >= ideal_mass_accreted[1]:
                mass_accreted_fractions.append(ideal_mass_accreted[1])
            else:
                mass_accreted_fractions.append(proj_layer_mass[1])
                mass_accreted_fractions[0] += (ideal_mass_accreted[1]-proj_layer_mass[1])
            
            compositions[targ_idx][2] = ((target_mass*target_core_frac)+(mass_accreted_fractions[0]))/largest_remnant_mass #new core frac for largest remnant
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2] #new mantle frac for largest remnant
            
            total_core_mass = (target_core_frac*target_mass)+(proj_core_frac*proj_mass) #how much total iron mass there is between the target and projectile
            largest_remnant_core_mass = largest_remnant_mass*compositions[targ_idx][2] #how much iron mass got into the largest remnant
            total_frag_core_mass = total_core_mass - largest_remnant_core_mass #how much iron mass is now left over in the fragments
            total_frag_mass = (target_mass+proj_mass) - largest_remnant_mass #total mass in fragments - THIS COULD CAUSE ERRORS MAYBE
            
            frag_core_frac = total_frag_core_mass/total_frag_mass #now gets the fraction of iron in each fragment by dividing the total fragment iron mass by the total fragment mass
            
            #If the remnant_core_frac is just barely above 1 or just barely below 0, it will be rounded down to 1 or up to 0
            if compositions[targ_idx][2] - 1 > 0 and compositions[targ_idx][2] - 1 < 1e-8: 
                compositions[targ_idx][2] = 1.0
                compositions[targ_idx][3] = 0.0
            if compositions[targ_idx][2] < 0 and compositions[targ_idx][2] > -1e-8:
                compositions[targ_idx][2] = 0.0
                compositions[targ_idx][3] = 1.0
            
            #If the frag_core_frac is just barely above 1 or just barely below 0, it will be rounded down to an even 1 or even 0
            if frag_core_frac - 1 > 0 and frag_core_frac - 1 < 1e-6: 
                frag_core_frac += 1.0 - frag_core_frac
            if frag_core_frac < 0 and frag_core_frac > -1e-6:
                frag_core_frac += 0.0 - frag_core_frac
            
            frag_mantle_frac = 1 - frag_core_frac
            
            #Error in case the fragments' core or mantle frac have a negative value
            if frag_core_frac < 0 or frag_mantle_frac < 0:
                print ('ERROR: Negative value for fragment composition encountered at', time)
                sys.exit(1)
                
                
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j], frag_core_frac, frag_mantle_frac] #creates a list filled with the necessary data for the fragments to go into the compositions array and the final output file - frags just given the composition of the projectile
                #Error in case a fragment has a negative mass
                if frag_masses[j] < 0:
                    print ('ERROR: Negative value for fragment mass encountered at', time)
                    sys.exit(1)
                compositions.append(frag_data)
            
            velocities.append(collision_velocities[i])
            if target_mass != largest_remnant_mass and sum(frag_masses) != proj_mass:
                disruptive_collision_velocities.append(float(collision_velocities[i]))
                disruptive_collision_target_masses.append(target_mass*obj_mass_conversion)
                disruptive_collision_proj_masses.append(proj_mass*obj_mass_conversion)
                disruptive_collision_lr_masses.append(largest_remnant_mass*obj_mass_conversion)
                disruptive_collision_target_cmfs.append(target_core_frac)
                disruptive_collision_fragment_masses.append([frag_masses[j]*obj_mass_conversion for j in range(len(frag_masses))])
                disruptive_collision_lr_cmfs.append(compositions[targ_idx][2])
            #FOLLOWING LISTS USED FOR FIGURE 3
                B_list.append(float(impact_parameters[i]))
                ideal_ejecta_CMF.append(ejecta_core_frac)
                if mass_accreted == 0:
                    actual_ejecta_CMF.append(0)
                else:
                    actual_ejecta_CMF.append(mass_accreted_fractions[0]/mass_accreted)
                collision_types.append('tab:blue')               
            
            #FOLLOWING LISTS USED FOR FIGURE 7
                if file_range_idx < 11:
                    efs.append('tab:blue')
                elif 11 <= file_range_idx < 21:
                    efs.append('tab:orange')
                elif 21 <= file_range_idx < 31:
                    efs.append('tab:green')
                elif 31 <= file_range_idx < 41:
                    efs.append('tab:purple')
                elif 41 <= file_range_idx < 51:
                    efs.append('tab:red')


################# PARTIAL EROSION & SUPER-CATASTROPHIC ##################      
                        
        if collision_type == 3 or collision_type == 4: #partial erosion, target abundances stay the same
            
            
            mass_lost = target_mass-largest_remnant_mass #change in mass of the target after the collision - this time mass is lost from target

            ejecta_core_frac = 0 #fraction of the eroded mass that is composed of core material (will depend on impact parameter)
            
            if target_core_radius == 0: #makes sure the core radius of the target isn't 0
                ejecta_core_frac = 0
            else:
                slope = (max_core_collision_frac-min_core_collision_frac)/(-2*target_core_radius)
                max_impact_parameter = proj_radius+target_core_radius #if impact parameter is less than this than it will strip off more core
                min_impact_parameter = proj_radius-target_core_radius #if impact parameter is greater than this than it will strip off less core
            #DETERMINES THE FRACTION OF CORE AND MANTLE THAT GETS LOST FROM TARGET
                if impact_parameter >= max_impact_parameter: #if the impact parameter is large (usually means oblique angle impact)
                    ejecta_core_frac += min_core_collision_frac #least amount of core is lost because cross-section of projectile won't intersect with core of target
                elif max_impact_parameter > impact_parameter > min_impact_parameter: #if it's between these extreme values
                    ejecta_core_frac += min_core_collision_frac+(slope*(impact_parameter-max_impact_parameter)) #it will strip off a certain amount of core according to this line equation
                elif impact_parameter <= min_impact_parameter: #if the impact parameter is small (often means a more head-on impact)
                    ejecta_core_frac += max_core_collision_frac #most amount of core is lost because cross-section of projectile will fully intersect with core of target
            
            
            ideal_mass_lost = [(mass_lost*ejecta_core_frac),(mass_lost*(1-ejecta_core_frac))] #first is core mass, second if mantle mass - how much of each layer would ideally be accreted if the projectile has the right composition
            target_layer_mass = [target_mass*target_core_frac, target_mass*(1-target_core_frac)] #how much mass is in each layer of projectile - core first, then mantle
            mass_lost_fractions = [] #will hold the actual fractions of each layer that get accreted from the projectile
            
            #these two statements make sure the correct amount of transferred mass is obtained from each layer - if there's not enough mass in a layer it takes the rest from the other layer
            if target_layer_mass[0] >= ideal_mass_lost[0]:
                mass_lost_fractions.append(ideal_mass_lost[0])
            else:
                mass_lost_fractions.append(target_layer_mass[0])
                ideal_mass_lost[1] += (ideal_mass_lost[0]-target_layer_mass[0])
                
            if target_layer_mass[1] >= ideal_mass_lost[1]:
                mass_lost_fractions.append(ideal_mass_lost[1])
            else:
                mass_lost_fractions.append(target_layer_mass[1])
                mass_lost_fractions[0] += (ideal_mass_lost[1]-target_layer_mass[1])
            
            
            compositions[targ_idx][2] = ((target_mass*target_core_frac)-(mass_lost_fractions[0]))/largest_remnant_mass #new core frac for largest remnant
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2] #new mantle frac for largest remnant
            
            total_core_mass = (target_core_frac*target_mass)+(proj_core_frac*proj_mass) #how much total iron mass there is between the target and projectile
            largest_remnant_core_mass = largest_remnant_mass*compositions[targ_idx][2] #how much iron mass got into the largest remnant
            total_frag_core_mass = total_core_mass - largest_remnant_core_mass #how much iron mass is now left over in the fragments
            total_frag_mass = (target_mass+proj_mass) - largest_remnant_mass #total mass in fragments - THIS COULD CAUSE ERRORS MAYBE
            
            frag_core_frac = total_frag_core_mass/total_frag_mass #now gets the fraction of iron in each fragment by dividing the total fragment iron mass by the total fragment mass
            
            #If the remnant_core_frac is just barely above 1 or just barely below 0, it will be rounded down to 1 or up to 0
            if compositions[targ_idx][2] - 1 > 0 and compositions[targ_idx][2] - 1 < 1e-8: 
                compositions[targ_idx][2] = 1.0
                compositions[targ_idx][3] = 0.0
            if compositions[targ_idx][2] < 0 and compositions[targ_idx][2] > -1e-8:
                compositions[targ_idx][2] = 0.0
                compositions[targ_idx][3] = 1.0
            
            #If the frag_core_frac is just barely above 1 or just barely below 0, it will be rounded down to 1 or up to 0
            if frag_core_frac - 1 > 0 and frag_core_frac - 1 < 1e-6: 
                frag_core_frac += 1.0 - frag_core_frac
            if frag_core_frac < 0 and frag_core_frac > -1e-6:
                frag_core_frac += 0.0 - frag_core_frac
            
            frag_mantle_frac = 1 - frag_core_frac
            
            #Error in case the fragments' core or mantle frac have a negative value
            if frag_core_frac < 0 or frag_mantle_frac < 0: 
                print(largest_remnant_mass)
                print(frag_core_frac)
                print ('ERROR: Negative value for fragment composition encountered at', time)
                sys.exit(1)
            
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j], frag_core_frac, frag_mantle_frac] #creates a list filled with the necessary data for the fragments to go into the compositions array and the final output file - frags just given the composition of the projectile
                #Error in case a fragment has a negative mass
                if frag_masses[j] < 0:
                    print ('ERROR: Negative value for fragment mass encountered at', time)
                    sys.exit(1)
                compositions.append(frag_data)
            
            velocities.append(collision_velocities[i])
            if target_mass != largest_remnant_mass and sum(frag_masses) != proj_mass:
                disruptive_collision_velocities.append(float(collision_velocities[i]))
                disruptive_collision_target_masses.append(target_mass*obj_mass_conversion)
                disruptive_collision_proj_masses.append(proj_mass*obj_mass_conversion)
                disruptive_collision_lr_masses.append(largest_remnant_mass*obj_mass_conversion)
                disruptive_collision_target_cmfs.append(target_core_frac)
                disruptive_collision_fragment_masses.append([frag_masses[j]*obj_mass_conversion for j in range(len(frag_masses))])
                disruptive_collision_lr_cmfs.append(compositions[targ_idx][2])
            #FOLLOWING LISTS USED FOR FIGURE 3
                B_list.append(float(impact_parameters[i]))
                ideal_ejecta_CMF.append(ejecta_core_frac)
                if mass_lost == 0:
                    actual_ejecta_CMF.append(0)
                else:  
                    actual_ejecta_CMF.append(mass_lost_fractions[0]/mass_lost)
                collision_types.append('tab:orange')
            
            #FOLLOWING LISTS USED FOR FIGURE 7
                if file_range_idx < 11:
                    efs.append('tab:blue')
                elif 11 <= file_range_idx < 21:
                    efs.append('tab:orange')
                elif 21 <= file_range_idx < 31:
                    efs.append('tab:green')
                elif 31 <= file_range_idx < 41:
                    efs.append('tab:purple')
                elif 41 <= file_range_idx < 51:
                    efs.append('tab:red')
                
        compositions[targ_idx][1] = largest_remnant_mass #changes mass of target to its post-collision mass
        
        #checks to make sure projectile isn't a second largest remnant before deleting it
        for hsh in frag_hashes:
            if hsh != proj_hash and hsh == frag_hashes[-1]:
                destroyed_object_hashes.append(proj_hash) #IMPORTANT CHANGE - removes projectile particle from the compositions list if it gets destroyed in collision  
            elif hsh == proj_hash:
                break
            else:
                continue
        
        if collision_type == 1: #Does the same thing as loop above - specifically for a merger 
            destroyed_object_hashes.append(proj_hash) #IMPORTANT CHANGE - removes projectile particle from the compositions list if it gets destroyed in collision
            
        #checks to see if any composition values in the target particle are negative
        for j in range(len(compositions[targ_idx])):
            if j == 1:
                if compositions[targ_idx][j] < 0: 
                    print ('ERROR: Negative value for target mass encountered at', time)
                    sys.exit(1)
            elif j > 1:
                if compositions[targ_idx][j] < 0:
                    print(compositions[targ_idx][j])
                    print ('ERROR: Negative value encountered in target data at', time)
                    sys.exit(1)
    
    #Destroys particles by hash      
    for hsh in destroyed_object_hashes:
        for j in range(len(compositions)):
            if compositions[j][0] == hsh:
                compositions.pop(j)
                break
            else:
                continue
    
    e_file = open(ejection_file, 'r')
    ejections_raw = [line.split() for line in e_file.readlines()]
    ejections = [int(ejections_raw[i][2]) for i in range(len(ejections_raw))]
    
    for hsh in ejections:
        for j in range(len(compositions)):
            if compositions[j][0] == hsh:
                compositions.pop(j)
                break
            else:
                continue
       
    e_file.close()

            
    return(compositions)

######## FILE WRITING FUNCTION ###########
#output has same format as input file
def write_output(compositions, composition_output_file):
    f = open(composition_output_file, "w")
    for body in compositions: #goes through each particle in compositions
        for i in range(len(body)): #iterates over each index in the individual particle list
                f.write(str(body[i]) + ' ')
        f.write('\n')#go to a new line and to a new particle
    f.close()
 
############# OUTPUT LOOP ####################
min_core_frac = 0.0 #mimimum fraction of core material in ejecta
max_core_frac = 1.0 #maximum fraction of core material in ejecta

collision_file_pw = "/Users/nferich/Desktop/anna_collision_reports/new_collision_reports/new_collision_report"
comp_input_file_pw = "/Users/nferich/Desktop/anna_collision_reports/mantle_stripping_input/uni_mantle_stripping_input"
impact_param_file_pw = "impact_parameters/impact_parameters"
ejec_file_pw = "ejections/ejections"
velocity_file_pw = "velocities/velocities"

file_range = np.arange(1,51,1)

start_time = time.time() #Timer to see how long running this code takes
############## MISCELLANEOUS GRAPHS #############

min_cmf = 0.0
max_cmf  = 1.0
for i in file_range:
    collision_file = collision_file_pw + str(i) + ".txt"
    comp_input_file = comp_input_file_pw + str(1) + ".txt"
    impact_param_file = impact_param_file_pw + str(i) + ".txt"
    ejec_file = ejec_file_pw + str(i) + ".txt"
    velo_file = velocity_file_pw + str(i) + ".txt"
    final_composition = track_composition(collision_file, comp_input_file, impact_param_file, ejec_file, velo_file, min_cmf, max_cmf, i)

"""for i in range(len(disruptive_collision_target_masses)):
    if disruptive_collision_target_masses[i] <= disruptive_collision_lr_masses[i]:
        collision_types[i] = 'tab:blue'
        print(i)
        print(disruptive_collision_target_masses[i])
        print(disruptive_collision_proj_masses[i])
        print(disruptive_collision_fragment_masses[i])"""

non_disruptive_collision_idxs = []
for i in range(len(disruptive_collision_target_masses)):
    if disruptive_collision_target_masses[i] == disruptive_collision_lr_masses[i]:
        if disruptive_collision_proj_masses[i] == sum(disruptive_collision_fragment_masses[i]):
            if len(disruptive_collision_fragment_masses[i]) == 1:
                non_disruptive_collision_idxs.append(i)
      
#if target_mass != largest_remnant_mass and sum(frag_masses) != proj_mass and len(frag_masses) == 1:
print(len(non_disruptive_collision_idxs))
print(len(disruptive_collision_lr_masses))
fig3, (ax3, ax4) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6,10))
fig3.subplots_adjust(hspace=0.03)
ax3.scatter(B_list, ideal_ejecta_CMF, color=collision_types, s=5.0, marker='o', linewidths=0, alpha=0.6)
ax3.set_ylabel('Ideal Ejecta CMF', fontsize='large')   
ax3.tick_params(length=3, width=1, which='major')
ax3.minorticks_on()
ax3_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='Accretive Collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='Erosive Collision', markerfacecolor='tab:orange', markersize=10.0)]
legend3 = ax3.legend(handles=ax3_legend_elements, loc = 'upper right', framealpha = .7)
ax4.scatter(B_list, actual_ejecta_CMF, color=collision_types, s=5.0, marker='o', linewidths=0, alpha=0.6)
ax4.set_xlabel('B/$R_{targ}$', fontsize='large')
ax4.set_ylabel('Actual Ejecta CMF', fontsize='large')   
ax4.minorticks_on()

plt.savefig('graphs/B_vs_ideal_and_actual_CMF.pdf', bbox_inches='tight', pad_inches=0.01)
    
fig7, ax7 = plt.subplots(figsize=(6,5))
ax7.scatter(B_list, ideal_ejecta_CMF, color=efs, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('B/$R_{targ}$', fontsize='large')
plt.ylabel('Ideal Ejecta CMF', fontsize='large')   
ax7.minorticks_on()
#plt.grid()
ax7_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='ef=3', markerfacecolor='tab:blue', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='ef=5', markerfacecolor='tab:orange', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:green', marker='o', lw=0.0, label='ef=7', markerfacecolor='tab:green', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:purple', marker='o', lw=0.0, label='ef=10', markerfacecolor='tab:purple', markersize=np.sqrt(.1*1000)),
                      Line2D([], [], color='tab:red', marker='o', lw=0.0, label='ef=15', markerfacecolor='tab:red', markersize=np.sqrt(.1*1000)),]
legend = plt.legend(handles=ax7_legend_elements, loc = 'upper right', framealpha = .7)
plt.gca().add_artist(legend)
plt.savefig('graphs/B_vs_ideal_CMF_with_ef.pdf', bbox_inches='tight', pad_inches=0.01)

fig8, ax8 = plt.subplots(figsize=(6,5))
ax8.scatter(B_list, disruptive_collision_velocities, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('B/$R_{targ}$', fontsize='large')
plt.ylabel('$v_{imp}/v_{esc}$', fontsize='large')
plt.yscale("log")  
plt.ylim([0.5,100])
ax8_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax8_legend_elements, loc = 'upper left', framealpha = .7)
plt.savefig('graphs/B_vs_vimp_over_v_esc.pdf', bbox_inches='tight', pad_inches=0.01)

fig9, ax9 = plt.subplots(figsize=(6,5))
ax9.scatter(disruptive_collision_target_masses, disruptive_collision_velocities, c=collision_types, s=1.5, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Target Mass ($M_{\u2295}$)', fontsize='large')
plt.ylabel('$v_{imp}/v_{esc}$', fontsize='large')
plt.yscale("log")  
plt.xscale("log")
ax9.axvline(.093, label = 'Embryo Mass', color = 'tab:red', linestyle = '--', alpha=.3)
ax9.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.3)
plt.ylim([0.5,100])
plt.legend(loc='upper right')
ax9_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax9_legend_elements, loc = 'upper right', framealpha = .7)
plt.savefig('graphs/target_mass_vs_vimp_over_vesc.pdf', bbox_inches='tight', pad_inches=0.01)

fig10, ax10 = plt.subplots(figsize=(6,5))
ax10.scatter(disruptive_collision_velocities, ideal_ejecta_CMF, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('$v_{imp}/v_{esc}$', fontsize='large')
plt.ylabel('Ideal Ejecta CMF', fontsize='large')
plt.xscale("log")
plt.xlim([0.5,100])
ax10_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax10_legend_elements, loc = 'upper left', framealpha = .7)
plt.savefig('graphs/vimp_over_vesc_vs_ideal_CMF.pdf', bbox_inches='tight', pad_inches=0.01)

fig11, ax11 = plt.subplots(figsize=(6,5))
ax11.scatter(disruptive_collision_target_masses, ideal_ejecta_CMF, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Target Mass ($M_{\u2295}$)', fontsize='large')
plt.ylabel('Ideal Ejecta CMF', fontsize='large')
plt.xscale("log")
ax11.axvline(.093, label = 'Embryo Mass', color = 'tab:red', linestyle = '--', alpha=.3)
ax11.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.3)
ax11_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax11_legend_elements, loc = 'upper left', framealpha = .7)
plt.savefig('graphs/targ_mass_vs_ideal_CMF.pdf', bbox_inches='tight', pad_inches=0.01)

total_radius = [calc_radius(disruptive_collision_target_masses[i]/obj_mass_conversion, 0, mantle_density, 0)+calc_radius(disruptive_collision_proj_masses[i]/obj_mass_conversion, 0, mantle_density, 0) for i in range(len(disruptive_collision_target_masses))] #supposed to have units of au
total_mass = [(disruptive_collision_target_masses[i]+disruptive_collision_proj_masses[i])/obj_mass_conversion for i in range(len(disruptive_collision_target_masses))]
v_esc_au_yr = [np.sqrt(2*G*total_mass[i]/total_radius[i]) for i in range(len(disruptive_collision_target_masses))]
v_esc = [v_esc_au_yr[i]*distance_conversion/(time_conversion*1000) for i in range(len(disruptive_collision_target_masses))] #should be km/s
v_impact = [disruptive_collision_velocities[i]*v_esc[i] for i in range(len(disruptive_collision_target_masses))]

fig12, (ax12, ax18) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6,10))
fig12.subplots_adjust(hspace=0.03)
ax12.scatter(v_impact, ideal_ejecta_CMF, c=collision_types, s=5.0, marker='o', linewidths=0, alpha=0.7)
ax12.set_ylabel('Ideal Ejecta CMF', fontsize='large')
plt.xscale("log")
ax12.minorticks_on()
ax12.set_xlim([1.0,400])
ax12_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend12 = ax12.legend(handles=ax12_legend_elements, loc = 'upper right', framealpha = .4)
ax18.scatter(v_impact, actual_ejecta_CMF, c=collision_types, s=5.0, marker='o', linewidths=0, alpha=0.7)
ax18.set_xlabel('$v_{imp}$ (km/s)', fontsize='large')
ax18.set_ylabel('Actual Ejecta CMF', fontsize='large')
plt.xscale("log")
ax18.set_xlim([1.0,400])
plt.savefig('graphs/vimp_vs_ideal_and_actual_ejecta_CMF.pdf', bbox_inches='tight', pad_inches=0.01)

fig13, ax13 = plt.subplots(figsize=(6,5))
ax13.scatter(disruptive_collision_target_masses, v_impact, c=collision_types, s=5.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Target Mass ($M_{\u2295}$)', fontsize='large')
plt.ylabel('$v_{imp}$ (km/s)', fontsize='large')
plt.yscale("log")  
plt.xscale("log")
ax13.axvline(.093, label='Embryo Mass', color='tab:red', linestyle='--', alpha=0.4)
ax13.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=0.4)
ax13.axvline(0.00465, label = 'Minimum Fragment Mass', color = 'tab:purple', linestyle = '--', alpha=0.4)
plt.ylim([0.5,250])
ax13_legend_elements = [Line2D([], [], color='tab:red', linestyle='--',label='Embryo Mass', alpha=0.4),
                      Line2D([], [], color='tab:green', linestyle='--', label='Planetesimal Mass', alpha=0.4),
                      Line2D([], [], color='tab:purple', linestyle='--', label='Minimum Fragment Mass', alpha=0.4)] 
ax13_legend_elements_2 = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='Accretive Collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='Erosive Collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax13_legend_elements, loc='lower right', framealpha=0.7)
legend2 = plt.legend(handles=ax13_legend_elements_2, loc = 'upper right', framealpha = 0.7)
plt.gca().add_artist(legend)
plt.gca().add_artist(legend2)
plt.savefig('graphs/target_mass_vs_v_imp.pdf', bbox_inches='tight', pad_inches=0.01)


fig14, ax14 = plt.subplots(figsize=(6,5))
ax14.scatter(disruptive_collision_target_masses, B_list, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Target Mass ($M_{\u2295}$)', fontsize='large')
plt.ylabel('B/$R_{targ}$', fontsize='large')
plt.xscale("log")
ax14.axvline(.093, label = 'Embryo Mass', color = 'tab:red', linestyle = '--', alpha=.3)
ax14.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.3)
ax14_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax14_legend_elements, loc = 'lower right', framealpha = .7)
plt.savefig('graphs/target_mass_vs_B.pdf', bbox_inches='tight', pad_inches=0.01)

fig15, ax15 = plt.subplots(figsize=(6,5))
ax15.scatter(disruptive_collision_lr_masses, v_impact, c=collision_types, s=1.5, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Largest Remnant Mass ($M_{\u2295}$)', fontsize='large')
plt.ylabel('$v_{imp}$ (km/s)', fontsize='large')
plt.yscale("log")  
plt.xscale("log")
ax15.axvline(.093, label = 'Embryo Mass', color = 'tab:red', linestyle = '--', alpha=.3)
ax15.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.3)
plt.ylim([0.5,1000])
ax15_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax15_legend_elements, loc = 'lower right', framealpha = .7)
plt.savefig('graphs/lr_mass_vs_vimp.pdf', bbox_inches='tight', pad_inches=0.01)

fig16, ax16 = plt.subplots(figsize=(6,5))
ax16.scatter(v_impact, disruptive_collision_lr_cmfs, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.7)
plt.xlabel('$v_{imp}$', fontsize='large')
plt.ylabel('Largest Remnant CMF', fontsize='large')
plt.xscale("log")
plt.xlim([1.0,1000])
ax16_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax16_legend_elements, loc = 'upper right', framealpha = .7)
plt.savefig('graphs/vimp_vs_lr_CMF.pdf', bbox_inches='tight', pad_inches=0.01)

fig17, ax17 = plt.subplots(figsize=(6,5))
ax17.scatter(disruptive_collision_target_masses, actual_ejecta_CMF, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Target Mass ($M_{\u2295}$)', fontsize='large')
plt.ylabel('Actual Ejecta CMF', fontsize='large')
plt.xscale("log")
ax17.axvline(.093, label = 'Embryo Mass', color = 'tab:red', linestyle = '--', alpha=.3)
ax17.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.3)
ax17_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax17_legend_elements, loc = 'lower right', framealpha = .7)
plt.savefig('graphs/target_mass_vs_actual_ejecta_cmf.pdf', bbox_inches='tight', pad_inches=0.01)

fig19, ax19 = plt.subplots(figsize=(6,5))
ax19.scatter(disruptive_collision_target_masses, disruptive_collision_lr_cmfs, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Target Mass ($M_{\u2295}$)', fontsize='large')
plt.ylabel('Largest Remnant CMF', fontsize='large')
plt.xscale("log")
ax19.axvline(.093, label = 'Embryo Mass', color = 'tab:red', linestyle = '--', alpha=.3)
ax19.axvline(.0093, label = 'Planetesimal Mass', color = 'tab:green', linestyle = '--', alpha=.3)
ax19_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax19_legend_elements, loc = 'upper right', framealpha = .7)
plt.savefig('graphs/target_mass_vs_lr_cmf.pdf', bbox_inches='tight', pad_inches=0.01)

fig20, ax20 = plt.subplots(figsize=(6,5))
ax20.scatter(disruptive_collision_target_cmfs, disruptive_collision_lr_cmfs, c=collision_types, s=3.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Target CMF', fontsize='large')
plt.ylabel('Largest Remnant CMF', fontsize='large')
ax20.minorticks_on()
ax20_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax20_legend_elements, loc = 'lower right', framealpha = .7)
plt.savefig('graphs/target_CMF_vs_LR_CMF.pdf', bbox_inches='tight', pad_inches=0.01)

fig21, ax21 = plt.subplots(figsize=(6,5))
ax21.scatter(B_list, v_impact, c=collision_types, s=2.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('B/$R_{targ}$', fontsize='large')
plt.ylabel('$v_{imp}$ (km/s)', fontsize='large')
plt.yscale("log")  
plt.ylim([1.0,1000])
ax8_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax8_legend_elements, loc = 'upper left', framealpha = .7)
plt.savefig('graphs/B_vs_vimp.pdf', bbox_inches='tight', pad_inches=0.01)

fig22, ax22 = plt.subplots(figsize=(6,5))
ax22.scatter(ideal_ejecta_CMF, actual_ejecta_CMF, c=collision_types, s=4.0, marker='o', linewidths=0, alpha=0.6)
plt.xlabel('Ideal Ejecta CMF', fontsize='large')  
plt.ylabel('Actual Ejecta CMF', fontsize='large') 
ax8_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax8_legend_elements, loc = 'upper left', framealpha = .7)
plt.savefig('graphs/ideal_CMF_vs_actual_CMF.pdf', bbox_inches='tight', pad_inches=0.01)

