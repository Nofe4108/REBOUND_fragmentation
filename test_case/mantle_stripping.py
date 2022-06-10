#!/usr/bin/env python

import numpy as np
import time
import sys

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
            if i >= 2:
                if line[i] != '&':
                    layer_comp.append(float(line[i]))
                if line[i] == '&':
                    particle_data.append(layer_comp)
                    layer_comp = []
                if i == len(line)-1:
                    particle_data.append(layer_comp)
                    layer_comp = []
        compositions.append(particle_data)
        
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
        
        if collision_type == 1: #perfect merger
            for i in range(no_layers): #index for each layer in target - (will change this once time dependence goes in)
                for j in range(len(last_target_abundances[i])): #index for how many elements are in a specific layer in target
                    compositions[targ_idx][i+2][j] = (float(last_target_abundances[i][j])*last_target_mass+float(last_projectile_abundances[i][j])*last_proj_mass)/target_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile mcomposition and mass

        if collision_type == 2: #partial accretion
            mass_accreted = target_mass-last_target_mass #change in mass of the target after the collision - this time mass is added to target
            print(mass_accreted)

track_composition()


print("--- %s seconds ---" % (time.time() - start_time))




###### EXTRA FUNCTIONS ######


