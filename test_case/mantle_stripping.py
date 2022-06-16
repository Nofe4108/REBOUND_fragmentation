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
        layer_comp = [] #will hold data about the composition of each layer in a particle before it's appended to the particle_data list
        particle_data = [] #This will hold data for individual particle - will contain its hash, mass, and lists containing the composition fractions for each layer
        for i in range(len(line)): #goes through the index for each element in a line
            if i == 0: #First value in each line is the hash so that'll automatically be added to the particle data list
                particle_data.append(int(line[i]))
            if i ==1: #Second value in each line is the particle mass so that'll automatically be added to the particle data list
                particle_data.append(float(line[i]))
            if i >= 2: #goes to the composition and layer parts of the particle
                if line[i] != '&': #if it's a composition value
                    layer_comp.append(float(line[i])) #add that to the layer list
                if line[i] == '&': #if it reaches a delimiter between layers
                    particle_data.append(layer_comp) #add the current layer to the data for the particle
                    layer_comp = [] #reset the list for the next layer
                if i == len(line)-1: #if it reaches the last composition value
                    particle_data.append(layer_comp) #add the final layer to the particle list
                    layer_comp = [] #reset the list for the next layer in next particle
        compositions.append(particle_data) #add all the data about the particle to the compositions list
        
    return(compositions)
    
    

start_time = time.time() #Timer to see how long running this code takes


###### MAIN FUNCTION #########
def track_composition(): #Main function that gives the compositions that will be outputted to the file
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/mantle_stripping_input.txt", 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    compositions = organize_compositions(init_compositions) #organzies the data from the file
    for line in compositions: #loops through each particle in the array
        if sum([sum(x) for x in line[2:]]) != 1.0: #makes sure composition fractions add up to one for each particle (uses second value in list and onwards)
            print ('ERROR: Realative abundances do not add up to 1.0')
            sys.exit(1)
                
    try: #init_hashes is the initial hashes in the sim as it goes through this loop? 
        init_hashes = [int(x[0]) for x in init_compositions] #makes sure hash input can be turned into integer, otherwise it assigns it a value?
    except:
        init_hashes = [x[0].value for x in init_compositions]
    
    no_species = sum(len(x) for x in compositions[0][2:]) #gives number of specie included in each layer - all lines in the input file must have the same number of elements so if an element in a certain body is 0 then it must be inputted as 0 and not left blank
    no_layers = len(compositions[0])-2#gives the number of layers that have been inputted - will remain constant
     
    file = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/collision_report.txt", 'r')
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
        last_target_abundances = compositions[targ_idx][-no_layers:] #list of composition fractions for target before the collision - [-no.species:] -no.species means it goes back from the end of the list to the versy first composition fraction and the colon means it adds all values from the list past that into a list 
        last_projectile_abundances = compositions[proj_idx][-no_layers:] #same thing as above except for projectile
        no_frags = int((len(block)-5)/2) #gives number of frags formed in collision - gets rid of 5 values from length which are basically about the collision type and target and projectile then divides in half since the list contains both the fragment hash and its mass
        frag_hashes = [int(block[i*2+3]) for i in range(1,no_frags+1)] #list of the hashes of the fragments - jumps from hash to hash using (i*2+3) - range length uses no_frags variable to get the correct amount of hashes
        frag_masses = [float(block[i*2+4]) for i in range(1,no_frags+1)]#same thing as above except for a list of the fragment masses 
        frag_mass = frag_masses[0]
          
                                               
 ######################## PERFECT MERGER ##########################
        if collision_type == 1: #perfect merger
            for i in range(no_layers): #index for each layer in target - (will change this once time dependence goes in)
                for j in range(len(last_target_abundances[i])): #index for how many elements are in a specific layer in target
                    compositions[targ_idx][i+2][j] = (float(last_target_abundances[i][j])*last_target_mass+float(last_projectile_abundances[i][j])*last_proj_mass)/target_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile mcomposition and mass
 
 ####################### PARTIAL ACCRETION ######################
   
        if collision_type == 2: #partial accretion
           
            mass_accreted = target_mass-last_target_mass #change in mass of the target after the collision - this time mass is added to target
            mass_loss_fractions = [.4, .6] #list inputted by user - tells what portion of the accreted mass comes from a certain layer of the projectile - first element is the innermost layer
            layer_mass_accreted = [mass_accreted*i for i in mass_loss_fractions] #mass from a specific layer accreted by the target from the projectile 
            
            last_target_layer_abundances = []
            last_projectile_layer_abundances = []
            
            #loop initializes the above two lists the correct dimensions
            for i in range(no_layers):
                zeroes = []
                for j in range(len(last_target_abundances[i])):
                    zeroes.append(0)
                last_target_layer_abundances.append(zeroes)
                last_projectile_layer_abundances.append(zeroes)
                zeroes = []
            
            # This loop will change the two lists above into lists that describe the fractional make up of each individual layer instead of the entire body overall
            for i in range(no_layers): 
                for j in range(len(last_target_abundances[i])):
                    last_target_layer_abundances[i][j] = last_target_abundances[i][j]/sum(last_target_abundances[i])
                    last_projectile_layer_abundances[i][j] = last_projectile_abundances[i][j]/sum(last_projectile_abundances[i])  

            for i in range(no_layers): #index for each layer in target - (will change this once time dependence goes in)
                for j in range(len(last_target_abundances[i])): #index for how many elements are in a specific layer in target
                    compositions[targ_idx][i+2][j]= ((float(last_target_abundances[i][j])*last_target_mass)+(last_projectile_layer_abundances[i][j]*layer_mass_accreted[i]))/target_mass #changes the compositional fraction of the element after the collision
    
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j]]+last_projectile_abundances #creates a list filled with the necessary data for the fragments to go into the compisitions array and the final output file - frags just given the composition of the projectile
                compositions.append(frag_data)
            
            

################# PARTIAL EROSION & SUPER-CATASTROPHIC ##################      
                        
        if collision_type == 3 or collision_type == 4: #partial errosion, target abundances stay the same
            
            mass_lost = last_target_mass-target_mass #change in mass of the target after the collision - this time mass is lost from target
            mass_loss_fractions = [.4, .6] #list inputted by user - tells what portion of the accreted mass comes from a certain layer of the projectile - first element is the innermost layer
            layer_mass_lost = [mass_lost*i for i in mass_loss_fractions] #mass from a specific layer eroded from the target
            last_target_layer_abundances = []
            last_projectile_layer_abundances = []
            
            #loop initializes the above two lists the correct dimensions
            for i in range(no_layers):
                zeroes = []
                for j in range(len(last_target_abundances[i])):
                    zeroes.append(0)
                last_target_layer_abundances.append(zeroes)
                last_projectile_layer_abundances.append(zeroes)
                zeroes = [] 
                
             # This loop will change the two lists above into lists that describe the fractional make up of each individual layer instead of the entire body overall
            for i in range(no_layers): 
                for j in range(len(last_target_abundances[i])):
                    last_target_layer_abundances[i][j] = last_target_abundances[i][j]/sum(last_target_abundances[i])
                    last_projectile_layer_abundances[i][j] = last_projectile_abundances[i][j]/sum(last_projectile_abundances[i])  
            
            frag_abundances = []
            
            #goes through and creates the post-collision abundances for the fragments
            for i in range(no_layers): 
                layer = []
                for j in range(len(last_target_abundances[i])):
                    layer.append(((last_target_layer_abundances[i][j]*layer_mass_lost[i])+(last_projectile_abundances[i][j]*last_proj_mass))/(mass_lost+last_proj_mass))
                frag_abundances.append(layer)
            
            #adds the necessary data about a fragment to a list and adds that 
            for i in range(no_frags):
                frag_data = [frag_hashes[i], frag_masses[i]]+frag_abundances
                compositions.append(frag_data)
                
            
        compositions[targ_idx][1] = target_mass #changes mass of target to its post-collision mass
        
        
        
        return(compositions)

######## FILE WRITING FUNCTION ###########
#output has same format as input file
def write_output(compositions):
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/mantle_stripping_output.txt", "w")
    for body in compositions: #goes through each particle in compositions
        for i in range(len(body)): #iterates over each index in the individual particle list
            if i <=1: #these two indices are the particle's hash and mass
                f.write(str(body[i]) + ' ')
            else: #now goes to the composition lists for each layer
                for comp in body[i]:#goes through each element in each layer
                    f.write(str(comp)+ ' ')
                if i != len(body)-1: #if this isn't the last layer
                    f.write('& ')#add the delimiter that goes between layers
                else:#when this is the last element
                    f.write('\n')#go to a new line and to a new particle
    f.close()
        
                    

track_composition()
write_output(track_composition())


print("--- %s seconds ---" % (time.time() - start_time))





