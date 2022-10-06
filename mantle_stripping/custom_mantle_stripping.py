#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:43:21 2022

@author: nferich
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 18:01:57 2022

@author: nferich
"""

#!/usr/bin/env python

import numpy as np
import math
import time
import sys

core_density = 7874.0 #kg m^-3 - density of iron
mantle_density = 3330.0 #kg m^3 -  approximate density of Earth's upper mantle

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

def calc_core_radius(mass, core_frac, core_density):
    core_mass = mass*core_frac
    core_radius = ((3*core_mass)/(4*core_density))**(1/3)
    return(core_radius)

def calc_radius(mass, core_frac, mantle_density, core_radius):
    mantle_frac = 1-core_frac
    mantle_mass = mass*mantle_frac
    radius = (core_radius**3+((3*mantle_mass)/(4*mantle_density)))**(1/3)
    return(radius)
      
#ratio = (2/np.pi)*math.atan((2-impact_param_over_Rtarg)*Q_over_Qstar/damping_factor)   
start_time = time.time() #Timer to see how long running this code takes


###### MAIN FUNCTION #########
def track_composition(): #Main function that gives the compositions that will be outputted to the file
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/disk_simulations/disk_mantle_stripping_input.txt", 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) #organzies the data from the file
    for line in compositions: #loops through each particle in the array
        if sum(line[2:]) != 1.0: #makes sure composition fractions add up to one for each particle (uses second value in list and onwards)
            print ('ERROR: Realative abundances do not add up to 1.0')
            sys.exit(1)
                
    try: 
        init_hashes = [int(x[0]) for x in init_compositions] 
    except:
        init_hashes = [x[0].value for x in init_compositions]
    
     
    file = open("/Users/nferich/GitHub/REBOUND_fragmentation/disk_simulations/collision_report.txt", 'r')
    blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
    blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list
    
    #START OF BIG LOOP
    for i in range(len(blocks)): #iterates through each value in blocks list - THIS IS A VERY BIG LOOP THAT CONTAINS ALL OF THE MASS TRANSFER DECISION MAKING
        block = blocks[i].split() #seperates each long string full of the collision data into its own list to be parsed through
        time = float(block[0]) #time of collision is first value
        collision_type = int(block[1]) #type of collision is second
        if collision_type == 0: #This is just an elastic bounce so it just continues to the next set of collision data in the blocks list
            continue
        
        impact_parameter = float(block[2]) #impact parameter of collision
        target_hash = int(block[3]) #hash of the target
        largest_remnant_mass = float(block[4]) #mass of the target after collision
        target_radius = float(block[5]) #radius of target before collision
        proj_hash = int(block[6]) #hash of projectile
        proj_radius = float(block[7]) #radius of projectile before collision
        targ_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==target_hash][0] #this searches through the compositions array to see where the target hash matches up with a hash in the array - variable saves this index so it can go right to the correct row for the target in the compositions array during future oprations
        proj_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==proj_hash][0]  #same thing as above except for projectile
        target_mass = float(compositions[targ_idx][1]) #gets the mass of target before collision - uses the targ_index variable to go to the right row
        proj_mass = float(compositions[proj_idx][1]) #same thing as above except for projectile
        target_core_frac = compositions[targ_idx][2] #list of composition fractions for target before the collision - [-no.species:] -no.species means it goes back from the end of the list to the versy first composition fraction and the colon means it adds all values from the list past that into a list 
        proj_core_frac = compositions[proj_idx][2] #same thing as above except for projectile
        no_frags = int((len(block)-8)/2) #gives number of frags formed in collision - gets rid of 5 values from length which are basically about the collision type and target and projectile then divides in half since the list contains both the fragment hash and its mass
        frag_hashes = [int(block[i*2+6]) for i in range(1,no_frags+1)] #list of the hashes of the fragments - jumps from hash to hash using (i*2+3) - range length uses no_frags variable to get the correct amount of hashes
        frag_masses = [float(block[i*2+7]) for i in range(1,no_frags+1)]#same thing as above except for a list of the fragment masses 
        
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
                
            compositions[targ_idx][2] = ((target_core_frac*target_mass)+(proj_core_frac*proj_mass))/largest_remnant_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile composition and mass
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
            
 
 ####################### PARTIAL ACCRETION ######################
   
        if collision_type == 2: #partial accretion
        
            mass_accreted = largest_remnant_mass-target_mass #change in mass of the target after the collision - this time mass is added to target
 
            min_core_collision_frac = 0.0 #minimum fraction of core material largest remnant gained in ejecta
            max_core_collision_frac = 0.30 #maximum fraction of core material largest remnant gained in ejecta
            ejecta_core_frac = 0 #fraction of the accreted mass that is composed of core material (will depend on impact parameter)
            
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

            #these two loops make sure the correct amount of transferred mass is obtained from each layer - if there's not enough mass in a layer it takes the rest from the other layer
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
                
            compositions[targ_idx][2] = ((target_mass*target_core_frac)+(mass_accreted_fractions[0]))/largest_remnant_mass
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]

            total_core_mass = (target_core_frac*target_mass)+(proj_core_frac*proj_mass) #how much total iron mass there is between the target and projectile
            largest_remnant_core_mass = largest_remnant_mass*compositions[targ_idx][2] #how much iron mass got into the largest remnant
            total_frag_core_mass = total_core_mass - largest_remnant_core_mass #how much iron mass is now left over in the fragments
            total_frag_mass = (target_mass+proj_mass) - largest_remnant_mass #total mass in fragments - THIS COULD CAUSE ERRORS MAYBE
            
            frag_core_frac = total_frag_core_mass/total_frag_mass #now gets the fraction of iron in each fragment by dividing the total fragment iron mass by the total fragment mass
            
            #If the frag_core_frac is just barely above 1 or just barely below 0, it will be rounded down to an even 1 or even 0
            if frag_core_frac - 1 > 0 and frag_core_frac - 1 < 1e-8: 
                frag_core_frac += 1.0 - frag_core_frac
            if frag_core_frac < 0 and frag_core_frac > -1e-8:
                frag_core_frac += 0.0 - frag_core_frac
            
            frag_mantle_frac = 1 - frag_core_frac
            
            #Error in case the fragments' core or mantle frac have a negative value
            if frag_core_frac < 0 or frag_mantle_frac < 0: 
                print ('ERROR: Negative value for fragment composition encountered at', time)
                sys.exit(1)
                
                
            for i in range(no_frags):
                frag_data = [frag_hashes[i], frag_masses[i], frag_core_frac, frag_mantle_frac] #creates a list filled with the necessary data for the fragments to go into the compositions array and the final output file - frags just given the composition of the projectile
                
                #Error in case a fragment has a negative mass
                if frag_masses[i] < 0:
                    print ('ERROR: Negative value for fragment mass encountered at', time)
                    sys.exit(1)
            
                compositions.append(frag_data)


################# PARTIAL EROSION & SUPER-CATASTROPHIC ##################      
                        
        if collision_type == 3 or collision_type == 4: #partial erosion, target abundances stay the same
            
            mass_lost = target_mass-largest_remnant_mass #change in mass of the target after the collision - this time mass is lost from target
            
            min_core_collision_frac = 0.0 #minimum fraction of core material largest remnant lost in ejecta
            max_core_collision_frac = 0.30 #maximum fraction of core material largest remnant lost in ejecta
            ejecta_core_frac = 0 #fraction of the eroded mass that is composed of core material (will depend on impact parameter)
            
            slope = (max_core_collision_frac-min_core_collision_frac)/(-2*target_core_radius)
            
            max_impact_parameter = proj_radius+target_core_radius #if impact parameter is less than this than it will strip off more core
            min_impact_parameter = proj_radius-target_core_radius #if impact parameter is greater than this than it will strip off less core
            
            #DETERMINES THE FRACTION OF CORE AND MANTLE THAT GETS LOST FROM TARGET
            if impact_parameter >= max_impact_parameter: #if the impact parameter is large (usually means oblique angle impact)
                ejecta_core_frac = min_core_collision_frac #least amount of core is lost because cross-section of projectile won't intersect with core of target
            elif max_impact_parameter > impact_parameter > min_impact_parameter: #if it's between these extreme values
                ejecta_core_frac = min_core_collision_frac+(slope*(impact_parameter-max_impact_parameter)) #it will strip off a certain amount of core according to this line equation
            elif impact_parameter <= min_impact_parameter: #if the impact parameter is small (often means a more head-on impact)
                ejecta_core_frac = max_core_collision_frac #most amount of core is lost because cross-section of projectile will fully intersect with core of target
            
            ideal_mass_lost = [(mass_lost*ejecta_core_frac),(mass_lost*(1-ejecta_core_frac))] #first is core mass, second if mantle mass - how much of each layer would ideally be accreted if the projectile has the right composition
            proj_layer_mass = [proj_mass*proj_core_frac, proj_mass*(1-proj_core_frac)] #how much mass is in each layer of projectile - core first, then mantle
            
            mass_lost_fractions = [] #will hold the actual fractions of each layer that get accreted from the projectile
            
            #these two loops make sure the correct amount of transferred mass is obtained from each layer - if there's not enough mass in a layer it takes the rest from the other layer
            if proj_layer_mass[0] >= ideal_mass_lost[0]:
                mass_lost_fractions.append(ideal_mass_lost[0])
            else:
                mass_lost_fractions.append(proj_layer_mass[0])
                ideal_mass_lost[1] += (ideal_mass_lost[0]-proj_layer_mass[0])
                
            if proj_layer_mass[1] >= ideal_mass_lost[1]:
                mass_lost_fractions.append(ideal_mass_lost[1])
            else:
                mass_lost_fractions.append(proj_layer_mass[1])
                mass_lost_fractions[0] += (ideal_mass_lost[1]-proj_layer_mass[1])
                
            
            compositions[targ_idx][2] = ((target_mass*target_core_frac)-(mass_lost_fractions[0]))/largest_remnant_mass
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
            
            total_core_mass = (target_core_frac*target_mass)+(proj_core_frac*proj_mass) #how much total iron mass there is between the target and projectile
            largest_remnant_core_mass = largest_remnant_mass*compositions[targ_idx][2] #how much iron mass got into the largest remnant
            total_frag_core_mass = total_core_mass - largest_remnant_core_mass #how much iron mass is now left over in the fragments
            total_frag_mass = (target_mass+proj_mass) - largest_remnant_mass #total mass in fragments - THIS COULD CAUSE ERRORS MAYBE
            
            frag_core_frac = total_frag_core_mass/total_frag_mass #now gets the fraction of iron in each fragment by dividing the total fragment iron mass by the total fragment mass
            
            #If the frag_core_frac is just barely above 1 or just barely below 0, it will be rounded down to an even 1 or even 0
            if frag_core_frac - 1 > 0 and frag_core_frac - 1 < 1e-8: 
                frag_core_frac += 1.0 - frag_core_frac
            if frag_core_frac < 0 and frag_core_frac > -1e-8:
                frag_core_frac += 0.0 - frag_core_frac
            
            frag_mantle_frac = 1 - frag_core_frac
            
            #Error in case the fragments' core or mantle frac have a negative value
            if frag_core_frac < 0 or frag_mantle_frac < 0: 
                print ('ERROR: Negative value for fragment composition encountered at', time)
                sys.exit(1)
            
            for i in range(no_frags):
                frag_data = [frag_hashes[i], frag_masses[i], frag_core_frac, frag_mantle_frac] #creates a list filled with the necessary data for the fragments to go into the compositions array and the final output file - frags just given the composition of the projectile
                #Error in case a fragment has a negative mass
                if frag_masses[i] < 0:
                    print ('ERROR: Negative value for fragment mass encountered at', time)
                    sys.exit(1)
                compositions.append(frag_data)
                
            
        compositions[targ_idx][1] = largest_remnant_mass #changes mass of target to its post-collision mass
        
        #If a projectile gets destroyed during a collision then this loop will remove it from the compositions list and not add it to the output file
        for hsh in frag_hashes:
            if hsh != proj_idx and hsh == frag_hashes[-1]:
                compositions.pop(proj_idx) #IMPORTANT CHANGE - removes projectile particle from the compositions list if it gets destroyed in collision
            elif hsh == proj_idx:
                break
            else:
                continue
        
        if collision_type == 1: #Does the same thing as loop above - specifically for a merger 
            compositions.pop(proj_idx) #IMPORTANT CHANGE - removes projectile particle from the compositions list if it gets destroyed in collision
            
        
        #checks to see if any composition values in the target particle are negative
        for i in range(len(compositions[targ_idx])):
            if i == 1:
                if compositions[targ_idx][i] < 0: 
                    print ('ERROR: Negative value for target mass encountered at', time)
                    sys.exit(1)
            elif i > 1:
                if compositions[targ_idx][i] < 0:
                    print(compositions[targ_idx][i])
                    print ('ERROR: Negative value encountered in target data at', time)
                    sys.exit(1)
            
    return(compositions)

######## FILE WRITING FUNCTION ###########
#output has same format as input file
def write_output(compositions):
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/mantle_stripping/custom_mantle_stripping_output.txt", "w")
    for body in compositions: #goes through each particle in compositions
        for i in range(len(body)): #iterates over each index in the individual particle list
                f.write(str(body[i]) + ' ')
        f.write('\n')#go to a new line and to a new particle
    f.close()
        

write_output(track_composition())

print("--- %s seconds ---" % (time.time() - start_time))
