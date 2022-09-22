#!/usr/bin/env python

import math
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
        particle_data = [] #This will hold data for individual particle - will contain its hash, mass, and lists containing the composition fractions for each layer
        for i in range(len(line)): #goes through the index for each element in a line
            if i == 0: #First value in each line is the hash so that'll automatically be added to the particle data list
                particle_data.append(int(line[i]))
            if i ==1: #Second value in each line is the particle mass so that'll automatically be added to the particle data list
                particle_data.append(float(line[i]))
                particle_data.append(0.0) #the third value in the compositions list will be the time at which the body formed
            if i == 2:
                particle_data.append(1) #Forth value describes how many layers the body has
            if i > 2: #goes to the composition and layer parts of the particle
                if line[i] != '|': #if it's a composition value and not the delimiter
                    particle_data.append(float(line[i])) #add that to the layer list
                if line[i] == '|': #if it reaches a delimiter between layers
                    particle_data[3] += 1 #add the current layer to the data for the particle
        compositions.append(particle_data) #add all the data about the particle to the compositions list
        
    return compositions
    

start_time = time.time() #Timer to see how long running this code takes




###### MAIN FUNCTION #########
def track_composition(): #Main function that gives the compositions that will be outputted to the file
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/mantle_stripping/mantle_stripping_input.txt", 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) #organzies the data from the file
    for line in compositions: #loops through each particle in the array
        if sum(line[4:]) != 1.0: #makes sure composition fractions add up to one for each particle (uses second value in list and onwards)
            print ('ERROR: Realative abundances do not add up to 1.0')
            sys.exit(1)
           
    try: #init_hashes is the initial hashes in the sim as it goes through this loop? 
        init_hashes = [int(x[0]) for x in init_compositions] #makes sure hash input can be turned into integer, otherwise it assigns it a value
    except:
        init_hashes = [x[0].value for x in init_compositions]
    
    
  ############### MUST BE ENTERED BY USER #####################
    no_species_in_layers = [1, 2, 1] #gives a list that says how many species there are in each layer - INPUTTED BY USER
    diff_timescale = 1e6 #how long it takes for a layer to be differentiated
  ##########################################################
    
    no_species = sum(no_species_in_layers) #gives number of specie included in each layer - all lines in the input file must have the same number of elements so if an element in a certain body is 0 then it must be inputted as 0 and not left blank
    max_no_layers = len(no_species_in_layers)#gives the number of layers that have been inputted - will remain constant
     
    file = open("/Users/nferich/GitHub/REBOUND_fragmentation/mantle_stripping/collision_report_3.txt", 'r')
    blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
    blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list
    for i in range(len(blocks)): #iterates through each value in blocks list - THIS IS A VERY BIG LOOP THAT CONTAINS ALL OF THE MASS TRANSFER DECISIONMAKING
        block = blocks[i].split() #seperates each long string full of the collision data into its own list to be parsed through
        time = float(block[0]) #time of collision is first value
        collision_type = int(block[1]) #type of collision is second
        if collision_type == 0: #This is just an elastic bounce so it just continues to the next set of collision data in the blocks list
            continue

        target_hash = int(block[2]) #hash of the target
        target_mass = float(block[3]) #mass of the target after collision
        projectile_hash = int(block[4]) #hash of projectile
        targ_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==target_hash][0] #this searches through the compositions array to see where the target hash matches up with a hash in the array - variable saves this index so it can go right to the correct row for the target in the compositions array during future oprations
        proj_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==projectile_hash][0]  #same thing as above except for projectile
        last_target_mass = float(compositions[targ_idx][1]) #gets the mass of target before collision - uses the targ_index variable to go to the right row
        last_proj_mass = float(compositions[proj_idx][1]) #same thing as above except for projectile
        last_target_abundances = compositions[targ_idx][-no_species:] #list of composition fractions for target before the collision - [-no.species:] -no.species means it goes back from the end of the list to the versy first composition fraction and the colon means it adds all values from the list past that into a list 
        last_proj_abundances = compositions[proj_idx][-no_species:] #same thing as above except for projectile
        no_frags = int((len(block)-5)/2) #gives number of frags formed in collision - gets rid of 5 values from length which are basically about the collision type and target and projectile then divides in half since the list contains both the fragment hash and its mass
        frag_hashes = [int(block[i*2+3]) for i in range(1,no_frags+1)] #list of the hashes of the fragments - jumps from hash to hash using (i*2+3) - range length uses no_frags variable to get the correct amount of hashes
        frag_masses = [float(block[i*2+4]) for i in range(1,no_frags+1)]#same thing as above except for a list of the fragment masses 
        targ_creation_time = compositions[targ_idx][2] #time in sim when target was created
        proj_creation_time = compositions[proj_idx][2] #time in sim when projectile was created
        targ_no_layers = compositions[targ_idx][3] #number of layers in target
        proj_no_layers = compositions[proj_idx][3] #number of layers in projectile
        
                                              
 ######################## PERFECT MERGER ##########################
         
        if collision_type == 1: #perfect merger
            for i in range(no_species): #index for each species in target 
                    compositions[targ_idx][i+4] = (float(last_target_abundances[i])*last_target_mass+float(last_proj_abundances[i])*last_proj_mass)/target_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile mcomposition and mass
 
 ####################### PARTIAL ACCRETION ######################
   
        if collision_type == 2: #partial accretion
           
            mass_accreted = target_mass-last_target_mass #change in mass of the target after the collision - this time mass is added to target
            
            time = 1.5e6
            
            #Following logic sees if the projectile has differentiated at all since its formation and the current collision it's a part of 
            if time - proj_creation_time > diff_timescale and proj_no_layers < max_no_layers: #if the projectile has been around longer than the differentiation timescale and isn't fully differentiated already
                proj_no_layers += math.floor((time - compositions[proj_idx][2])/diff_timescale) #add a number layers proportional to the floored integer of the age of the projectile divided by the differentiation timescale 
                if proj_no_layers > max_no_layers: #if the number of layers in the projectile is now over the max number
                    proj_no_layers = max_no_layers #change the number of layers to the max number
            
            layered_last_proj_abundances = [] #list that will contained the abundances of the projectile split into layers     
            
            #Following logic will split the abundances into a nested list where the sublists represent layers that the projectile has
            layer = [] #holds a single layer before appending
            comp_counter = 0 #this will keep track of the proper starting index for a new layer
            for i in range(proj_no_layers): #for number of layers the projectile has
                for j in range(no_species_in_layers[i]): #for the number of species that are supposed to be in a particular layer i 
                    layer.append(last_proj_abundances[comp_counter+j]) #adds the proper element to the layer 
                layered_last_proj_abundances.append(layer) #adds the completed layer to the layered abundance list
                layer = [] #resets layer list for new layer
                comp_counter += no_species_in_layers[i] #moves the starting index up to the next element in the abundances list
            if proj_no_layers < max_no_layers: #if the number of layers in the projectile is less than the max numbers
                remaining_composition = last_proj_abundances[comp_counter:] #creates a list filled with the remaining abundances that will just be placed in the top layer
                for comp in remaining_composition: 
                    layered_last_proj_abundances[-1].append(comp) #adds the remaining leftover elements to the topmost (which in the list is considered the rightmost) layer
            
            proj_layer_mass_fractions = [sum(i) for i in layered_last_proj_abundances] #percentage of mass contained in each layer of projectile
            
            mass_accreted_single_frac = mass_accreted*((proj_no_layers**2)+proj_no_layers)/2 #Addition factorial of the number of layers currently in the projectile - will be used to determine mass lost from each layer - if 2 layers ratio will be 1 to 2 - if 3 layers ratio will be 1 to 2 to 3 - etc.
            ideal_layer_mass_accreted = [mass_accreted_single_frac*(i+1) for i in range(proj_no_layers)] #the ideal amount of mass that would come from each layer in a collision - may not match up correctly
            proj_layer_mass_accreted = [0 for i in range(proj_no_layers)]
            
            for i in range(proj_no_layers):
                if ideal_layer_mass_accreted[i] > proj_layer_mass_fractions[i] and : #if there's not enough mass to take out of the layer
                    leftover_mass = ideal_layer_mass_accreted[i] - proj_layer_mass_fractions[i]
                elif
                
                else:
                    proj_layer_mass_accreted[i] = ideal_layer_mass_accreted[i]
                    
                    
                    
            
            proj_layer_mass_accreted = []
            for i in range(proj_no_layers):
                if i == proj_no_layers-1:
                    for i in range()
                
            
            proj_layer_mass_accreted = [mass_accreted*i for i in proj_layer_mass_fractions] #mass from a specific layer accreted by the target from the projectile - directly proportional to the amount of mass contained in each layer

            proj_relative_layer_abundances = []
            
            #loop initializes the list above to the correct dimensions
            for i in range(compositions[proj_idx][3]):
                zeroes = []
                for j in range(len(layered_last_proj_abundances[i])):
                    zeroes.append(0.0)
                proj_relative_layer_abundances.append(zeroes)
                zeroes = []
            
            # This logic will change the the list above into a list that describes the fractional make up of each individual layer instead of the entire body overall
            for i in range(proj_no_layers): 
                for j in range(len(layered_last_proj_abundances[i])):
                    if sum(layered_last_proj_abundances[i]) == 0: #if there's no mass in a layer
                        proj_relative_layer_abundances[i][j] = 0 #then the relative abundances of that layer are just 0
                    else:
                        proj_relative_layer_abundances[i][j] = layered_last_proj_abundances[i][j]/sum(layered_last_proj_abundances[i])  
            
            species_mass_accreted = [] #list will hold how much mass from each indicidual species was accreted by the target
            
            #This logic takes the list containing how much mass was accreted from each layer of the projectile and the list containing the abundances of each individual layer to figure out how much mass from each species the target gained
            for i in range(proj_no_layers): #for number of layers in projectile
                for j in range(len(proj_relative_layer_abundances[i])): #for number of elements in each layer
                    species_mass_accreted.append(proj_layer_mass_accreted[i]*proj_relative_layer_abundances[i][j]) #add the mass accreted for that particular species to the list
            
            for i in range(no_species):
                compositions[targ_idx][i+4] = ((float(last_target_abundances[i])*last_target_mass)+(species_mass_accreted[i]))/target_mass #changes the compositional fraction of the element after the collision

            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j], time, 1]+last_proj_abundances #creates a list filled with the necessary data for the fragments to go into the compisitions array and the final output file - frags just given the composition of the projectile - assumes fragments are undifferentiated
                compositions.append(frag_data)
            
            #I don't think an error is necessary for this type of collision - covered by the error statement at the very end of the function
           

################# PARTIAL EROSION & SUPER-CATASTROPHIC ##################      
                        
        if collision_type == 3 or collision_type == 4: #partial errosion, target abundances stay the same
            
            mass_lost = last_target_mass-target_mass #change in mass of the target after the collision - this time mass is lost from target

            
            #Following logic sees if the target has differentiated at all since its formation and the current collision it's a part of 
            if time - targ_creation_time > diff_timescale and targ_no_layers < max_no_layers: #if the target has been around longer than the differentiation timescale and isn't fully differentiated already
                targ_no_layers += math.floor((time - compositions[targ_idx][2])/diff_timescale) #add a number layers proportional to the floored integer of the age of the target divided by the differentiation timescale 
                if targ_no_layers > max_no_layers: #if the number of layers in the target is now over the max number
                    targ_no_layers = max_no_layers #change the number of layers to the max number
            
            layered_last_targ_abundances = []
            
            #Following logic will split the abundances into a nested list where the sublists represent layers that the target has
            layer = [] #holds a single layer before appending
            comp_counter = 0 #this will keep track of the proper starting index for a new layer
            for i in range(targ_no_layers): #for number of layers the projectile has
                for j in range(no_species_in_layers[i]): #for the number of species that are supposed to be in a particular layer i 
                    layer.append(last_target_abundances[comp_counter+j]) #adds the proper element to the layer 
                layered_last_targ_abundances.append(layer) #adds the completed layer to the layered abundance list
                layer = [] #resets layer list for new layer
                comp_counter += no_species_in_layers[i] #moves the starting index up to the next element in the abundances list
            if targ_no_layers < max_no_layers: #if the number of layers in the target is less than the max numbers
                remaining_composition = last_target_abundances[comp_counter:] #creates a list filled with the remaining abundances that will just be placed in the top layer
                for comp in remaining_composition: 
                    layered_last_targ_abundances[-1].append(comp) #adds the remaining leftover elements to the topmost (which in the list is considered the rightmost) layer
                    
            targ_layer_mass_fractions = [sum(i) for i in layered_last_targ_abundances] #percentage of mass contained in each layer of projectile
            targ_layer_mass_lost = [mass_lost*i for i in targ_layer_mass_fractions] #mass from a specific layer accreted by the target from the projectile - directly proportional to the amount of mass contained in each layer

            targ_relative_layer_abundances = []
            
            #loop initializes the list above to the correct dimensions
            for i in range(compositions[targ_idx][3]):
                zeroes = []
                for j in range(len(layered_last_targ_abundances[i])):
                    zeroes.append(0.0)
                targ_relative_layer_abundances.append(zeroes)
                zeroes = []
                
            # This logic will change the the list above into a list that describes the fractional make up of each individual layer instead of the entire body overall
            for i in range(targ_no_layers): 
                for j in range(len(layered_last_targ_abundances[i])):
                    if sum(layered_last_targ_abundances[i]) == 0: #if there's no mass in a layer
                        targ_relative_layer_abundances[i][j] = 0 #then the relative abundances of that layer are just 0
                    else:
                        targ_relative_layer_abundances[i][j] = layered_last_targ_abundances[i][j]/sum(layered_last_targ_abundances[i]) 
                    
                    
            frag_abundances = []
            species_mass_lost = [] #list will hold how much mass from each indicidual species was lost from the target
            
            #This logic takes the list containing how much mass was lost from each layer of the target and the list containing the abundances of each individual layer to figure out how much mass from each species the target lost
            for i in range(targ_no_layers): #for number of layers in projectile
                for j in range(len(targ_relative_layer_abundances[i])): #for number of elements in each layer
                    species_mass_lost.append(targ_layer_mass_lost[i]*targ_relative_layer_abundances[i][j]) #add the mass accreted for that particular species to the list
            print(species_mass_lost)
            print(last_proj_abundances)
            print(last_proj_mass)
            print(no_species)
            
            #goes through and creates the post-collision abundances for the fragments
            for i in range(no_species): 
                frag_abundances.append((species_mass_lost[i]+(last_proj_abundances[i]*last_proj_mass))/(mass_lost+last_proj_mass))
            
            #checks to see if there are any negative values in the frag data - only care about the abundances,hash and mass of the fragments should be fine
            for species in frag_abundances:
                if any(float(i) < 0 for i in layer) == True:
                    print ('ERROR: Negative value encountered in frag data at', time)
                    sys.exit(1)
            
            #adds the necessary data about a fragment to a list and adds that to the main compositions list as a new particle
            for i in range(no_frags):
                frag_data = [frag_hashes[i], frag_masses[i], time, 1]+frag_abundances # assumes fragments are undifferentiated
                compositions.append(frag_data)
                      
        
        compositions[targ_idx][1] = target_mass #changes mass of target to its post-collision mass
        
        #checks to see if any composition values in the target particle are negative
        for i in range(len(compositions[targ_idx])):
            if i == 1:
                if compositions[targ_idx][i] < 0: 
                    print ('ERROR: Negative value for target mass encountered at', time)
                    sys.exit(1)
            elif i > 3:
                if compositions[targ_idx][i] < 0:
                    print ('ERROR: Negative value encountered in target data at', time)
                    sys.exit(1)
            
    return(compositions)

######## FILE WRITING FUNCTION ###########
#output has same format as input file
def write_output(compositions):
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/mantle_stripping/better_mantle_stripping_output.txt", "w")
    for body in compositions: #goes through each particle in compositions
        for i in range(len(body)): #iterates over each index in the individual particle list
            if i <=3: #these two indices are the particle's hash and mass
                f.write(str(body[i]) + ' ')
            else: #now goes to the composition for each layer
                f.write(str(body[i])+ ' ')
                if i == len(body)-1: #if this isn't the last layer
                    f.write('\n')#go to a new line and to a new particle
    f.close()
        

write_output(track_composition())

print("--- %s seconds ---" % (time.time() - start_time))
