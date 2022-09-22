#!/usr/bin/env python

import numpy as np
import time
import sys


######## COMPOSITION DATA ORGANIZING FUNCTION ###########
#Function takes the data from the composition input file and organizes it in a large array
#The created list has seperate values that reresent the particles
#within these nested lists will be the hash and mass of the particles
#This will then be followed by more nested lists which represent the layers of the particle with its fractional compositions
def organize_compositions(init_compositions):
    compositions = []
    for line in init_compositions: #iterates through each line from file
        particle_data = [] #This will hold data for individual particle - will contain its hash, mass, and lists containing the fractions for each layer
        for i in range(len(line)): #goes through the index for each element in a line
            if i == 0: #First value in each line is the hash so that'll be an integer
                particle_data.append(int(line[i]))
            else: #The other values will be floats
                particle_data.append(float(line[i]))
        compositions.append(particle_data) #add all the data about the particle to the compositions list
    return(compositions)

#Variables needed for the iron frac function
e =2.71828

aFe = 0.2099
bFe = 2.1516
bFeSc = 6.4890

qSc = 1.26 #This is the value of Qr/Qrd in Reinhardt et al. (2022) that tells you the cutoff for a supercastastrophic collision
mass_ratio_cutoff = 1-(qSc/2) #This is the previous cutoff but now it terms of Mlr/Mtot

#BASICALLY EQ. 6 FROM REINHARDT ET AL. 2022
def iron_frac_equation(mass_lr, target_mass, proj_mass, target_fe_frac, proj_fe_frac):
    
    total_mass = target_mass+proj_mass #total mass involved in collision
    total_iron_frac = ((target_fe_frac*target_mass)+(proj_fe_frac*proj_mass))/(total_mass) #fraction of total mass that's iron
    energy_ratio = 2*(1-(mass_lr/total_mass)) #This is equivalent to Qr/Qrd in Reinhardt et al. (2022)
    
    if mass_lr/(target_mass+proj_mass) > mass_ratio_cutoff: #if it's not a supercatastrophic collision
        
        largest_remnant_iron_frac = total_iron_frac + aFe*(energy_ratio)**bFe
        return largest_remnant_iron_frac
    
    else:

        aFeSc = 1 - (total_iron_frac + aFe*(qSc**bFe))
        largest_remnant_iron_frac = 1 - aFeSc*e**(-bFeSc*(energy_ratio-qSc))
        return largest_remnant_iron_frac
    

start_time = time.time() #Timer to see how long running this code takes


###### MAIN FUNCTION #########
def track_composition(): #Main function that gives the compositions that will be outputted to the file
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/iron_frac_input.txt", 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) #organzies the data from the file
    for line in compositions: #loops through each particle in the array
        if sum(line[2:]) != 1.0: #makes sure composition fractions add up to one for each particle (uses second value in list and onwards)
            print ('ERROR: Realative abundances do not add up to 1.0')
            sys.exit(1)
             
    
    try: #init_hashes is the initial hashes in the sim as it goes through this loop? 
        init_hashes = [int(x[0]) for x in init_compositions] #makes sure hash input can be turned into integer, otherwise it assigns it a value
    except:
        init_hashes = [x[0].value for x in init_compositions]
     
    file = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/collision_report.txt", 'r')
    blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
    blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list
    
    #START OF BIG LOOP
    for i in range(len(blocks)): #iterates through each value in blocks list - THIS IS A VERY BIG LOOP THAT CONTAINS ALL OF THE MASS TRANSFER DECISIONMAKING
        block = blocks[i].split() #seperates each long string full of the collision data into its own list to be parsed through
        time = float(block[0]) #time of collision is first value
        collision_type = int(block[1]) #type of collision is second
        if collision_type == 0: #This is just an elastic bounce so it just continues to the next set of collision data in the blocks list
            continue

        target_hash = int(block[2]) #hash of the target
        largest_remnant_mass = float(block[3]) #mass of the target after collision
        projectile_hash = int(block[4]) #hash of projectile
        targ_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==target_hash][0] #this searches through the compositions array to see where the target hash matches up with a hash in the array - variable saves this index so it can go right to the correct row for the target in the compositions array during future oprations
        proj_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==projectile_hash][0]  #same thing as above except for projectile
        last_target_mass = float(compositions[targ_idx][1]) #gets the mass of target before collision - uses the targ_index variable to go to the right row
        last_proj_mass = float(compositions[proj_idx][1]) #same thing as above except for projectile
        last_target_iron_frac = compositions[targ_idx][2] #list of composition fractions for target before the collision - [-no.species:] -no.species means it goes back from the end of the list to the versy first composition fraction and the colon means it adds all values from the list past that into a list 
        last_proj_iron_frac = compositions[proj_idx][2] #same thing as above except for projectile
        no_frags = int((len(block)-5)/2) #gives number of frags formed in collision - gets rid of 5 values from length which are basically about the collision type and target and projectile then divides in half since the list contains both the fragment hash and its mass
        frag_hashes = [int(block[i*2+3]) for i in range(1,no_frags+1)] #list of the hashes of the fragments - jumps from hash to hash using (i*2+3) - range length uses no_frags variable to get the correct amount of hashes
        frag_masses = [float(block[i*2+4]) for i in range(1,no_frags+1)]#same thing as above except for a list of the fragment masses 

          
                                               
 ######################## PERFECT MERGER ##########################
        if collision_type == 1: #perfect merger
            
            compositions[targ_idx][2] = ((last_target_iron_frac*last_target_mass)+(last_proj_iron_frac*last_proj_mass))/largest_remnant_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile mcomposition and mass
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
 
 ####################### PARTIAL ACCRETION ######################
   
        if collision_type == 2: #partial accretion
            
            compositions[targ_idx][2] = iron_frac_equation(largest_remnant_mass, last_target_mass, last_proj_mass, last_target_iron_frac, last_proj_iron_frac)
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
            
            total_iron_mass = (last_target_iron_frac*last_target_mass)+(last_proj_iron_frac*last_proj_mass) #how much total iron mass there is between the target and projectile
            largest_remnant_iron_mass = largest_remnant_mass*compositions[targ_idx][2] #how much iron mass got into the largest remnant
            total_frag_iron_mass = total_iron_mass - largest_remnant_iron_mass #how much iron mass is now left over in the fragments
            
            frag_iron_frac = total_frag_iron_mass/sum(frag_masses) #now gets the fraction of iron in each fragment by dividing the total fragment iron mass by the total fragment mass
            frag_mantle_frac = 1 - frag_iron_frac
    
            for i in range(no_frags):
                frag_data = [frag_hashes[i], frag_masses[i], frag_iron_frac, frag_mantle_frac] #creates a list filled with the necessary data for the fragments to go into the compositions array and the final output file - frags just given the composition of the projectile
                compositions.append(frag_data)
        
            """mass_accreted = largest_remnant_mass-last_target_mass
            
            compositions[targ_idx][2] = ((last_target_mass*last_target_iron_frac)+(mass_accreted*last_proj_iron_frac))/largest_remnant_mass #iron frac of target post collision - the mass accreted from projectile has same iron frac as the projectile  
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
            
            total_iron_mass = (last_target_iron_frac*last_target_mass)+(last_proj_iron_frac*last_proj_mass) #how much total iron mass there is between the target and projectile
            largest_remnant_iron_mass = largest_remnant_mass*compositions[targ_idx][2] #how much iron mass got into the largest remnant
            total_frag_iron_mass = total_iron_mass - largest_remnant_iron_mass #how much iron mass is now left over in the fragments
            
            frag_iron_frac = total_frag_iron_mass/sum(frag_masses) #now gets the fraction of iron in each fragment by dividing the total fragment iron mass by the total fragment mass
            frag_mantle_frac = 1 - frag_iron_frac
    
            for i in range(no_frags):
                frag_data = [frag_hashes[i], frag_masses[i], frag_iron_frac, frag_mantle_frac] #creates a list filled with the necessary data for the fragments to go into the compisitions array and the final output file - frags just given the composition of the projectile
                compositions.append(frag_data)"""
                
            #I don't think an error is necessary for this type of collision - covered by the error statement at the very end of the function
           
################# PARTIAL EROSION & SUPER-CATASTROPHIC ##################      
                        
        if collision_type == 3 or collision_type == 4: #partial errosion, target abundances stay the same
        
            compositions[targ_idx][2] = iron_frac_equation(largest_remnant_mass, last_target_mass, last_proj_mass, last_target_iron_frac, last_proj_iron_frac)
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
            
            total_iron_mass = (last_target_iron_frac*last_target_mass)+(last_proj_iron_frac*last_proj_mass) #how much total iron mass there is between the target and projectile
            largest_remnant_iron_mass = largest_remnant_mass*compositions[targ_idx][2] #how much iron mass got into the largest remnant
            total_frag_iron_mass = total_iron_mass - largest_remnant_iron_mass #how much iron mass is now left over in the fragments
            
            frag_iron_frac = total_frag_iron_mass/sum(frag_masses) #now gets the fraction of iron in each fragment by dividing the total fragment iron mass by the total fragment mass
            frag_mantle_frac = 1 - frag_iron_frac
    
            for i in range(no_frags):
                frag_data = [frag_hashes[i], frag_masses[i], frag_iron_frac, frag_mantle_frac] #creates a list filled with the necessary data for the fragments to go into the compisitions array and the final output file - frags just given the composition of the projectile
                compositions.append(frag_data)
        
        
        compositions[targ_idx][1] = largest_remnant_mass #changes mass of target to its post-collision mass
        compositions.pop(proj_idx) #IMPORTANT CHANGE - removes projectile particle from the compositions list
        
        print(compositions)
        #checks to see if any composition values in the target particle are negative
        for i in range(len(compositions[targ_idx])):
            if i == 1:
                if compositions[targ_idx][i] < 0: 
                    print ('ERROR: Negative value for target mass encountered at', time)
                    sys.exit(1)
            elif i > 1:
                if compositions[targ_idx][i] < 0:
                    print ('ERROR: Negative value encountered in target data at', time)
                    sys.exit(1)
            
    return(compositions)

######## FILE WRITING FUNCTION ###########
#output has same format as input file
def write_output(compositions):
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/iron_frac_output.txt", "w")
    for body in compositions: #goes through each particle in compositions
        for i in range(len(body)): #iterates over each index in the individual particle list
                f.write(str(body[i]) + ' ')
        f.write('\n')#go to a new line and to a new particle
    f.close()
        

write_output(track_composition())

print("--- %s seconds ---" % (time.time() - start_time))





