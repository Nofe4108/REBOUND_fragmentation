#!/usr/bin/env python

import numpy as np
import time
import sys

start_time = time.time() #Timer to see how long running this code takes
def track_composition(): #Main function that gives the compositions that will be outputted to the file
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/test_composition_input.txt", 'r')
    init_compositions = [line.split() for line in f.readlines()] #reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
    for line in init_compositions: #loops through each list in the array 
        try:
            sum([float(x) for x in line[2:]]) == 1.0 #makes sure composition fractions add up to one for each particle (uses second value in list and onwards)
        except:
            print ('ERROR: Realative abundances do not add up to 1.0')
            sys.exit(1)
    try: #init_hashes is the initial hashes in the sim as it goes through this loop? 
        init_hashes = [int(x[0]) for x in init_compositions] #makes sure hash input can be turned into integer, otherwise it assigns it a value?
    except:
        init_hashes = [x[0].value for x in init_compositions]
    no_species = len(init_compositions[0])-2 #gives number of specie included in the list - all lines in the list must have the same length so if the calue of an element in a certain body is 0 then it must be inputted as 0 and not left blank (I CAN DO THE SAME THING USING A TAB TO INDICATE A NEW LAYER)
    file = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/collision_report.txt", 'r')
    blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
    blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list
    compositions = init_compositions #name change for initial compositions
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
        last_projectile_abundances = compositions[proj_idx][-no_species:] #same thing as above except for projectile
        no_frags = int((len(block)-5)/2) #gives number of frags formed in collision - gets rid of 5 values from length which are basically about the collision type and target and projectile then divides in half since the list contains both the fragment hash and its mass
        frag_hashes = [int(block[i*2+3]) for i in range(1,no_frags+1)] #list of the hashes of the fragments - jumps from hash to hash using (i*2+3) - range length uses no_frags variable to get the correct amount of hashes
        frag_masses = [float(block[i*2+4]) for i in range(1,no_frags+1)]#same thing as above except for a list of the fragment masses
        if collision_type == 1: #perfect merger
            for i in range(no_species):
                compositions[targ_idx][i+2]=(float(last_target_abundances[i])*last_target_mass+float(last_projectile_abundances[i])*last_proj_mass)/target_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile mcomposition and mass
        if collision_type == 2: #partial accretion
            mass_accreted = target_mass-last_target_mass #change in mass of the target after the collision - this time mass is added to target
            for i in range(no_species):
                compositions[targ_idx][i+2]=(float(last_target_abundances[i])*last_target_mass+float(last_projectile_abundances[i])*mass_accreted)/target_mass #changes the composition fraction of each specie in target - again a type of weighted average with the initial comp and mass of the target and the comp of the projectile and mass accreted
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j]]+last_projectile_abundances #creates a list filled with the necessary data for the fragments to go into the compisitions array and the final output file - frags just given the composition of the projectile
                compositions.append(frag_data) #adds the new frag data to the compositions array
                for i in range(len(frag_data)):  #THIS LOOP TESTS FOR NEGATIVE VALUES IN FRAG DATA - TRY/EXCEPT LOOP WASN"T WORKING FOR ME
                    if float(frag_data[i]) < 0:
                        print ('ERROR: Negative value encountered in frag data at', time)
                        sys.exit(1)
                """try:
                     any(n < 0 for n in frag_data) == False
                except:
                    print ('ERROR: Negative value encountered in frag data at', time)
                    sys.exit(1)"""
        if collision_type == 3 or collision_type == 4: #partial errosion, target abundances stay the same
            mass_lost = last_target_mass-target_mass #change in mass of the target after the collision - this time mass is lost from target
            frag_abundances = [(float(last_target_abundances[i])*mass_lost+float(last_projectile_abundances[i])*last_proj_mass)/np.sum(frag_masses) for i in range(no_species)] #determines abundances for fragments and puts it in a list with a value for each fragment - gives same value to all fragments - again a type of weighted value between the comp of the target and the mass lost from it and the comp and mass of the projectile 
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j]]+frag_abundances #creates a list filled with the necessary data for the fragments to go into the compisitions array and the final output file 
                try:
                     any(n < 0 for n in frag_data) == False
                except:
                    print ('ERROR: Negative value encountered in frag data at', time)
                    sys.exit(1)
                compositions.append(frag_data) #adds fragment data to the compositions array
        compositions[targ_idx][1]=target_mass #changes the mass of the target in the composition array to the post-collision value
        for i in range(len(compositions[targ_idx])):  #THIS LOOP TESTS FOR NEGATIVE VALUES IN COMPOSITIONS
                    if float(compositions[targ_idx][i]) < 0:
                        print ('ERROR: Negative value encountered at', time)
                        sys.exit(1)
        """try: 
            any(n < 0 for n in compositions[targ_idx]) == False
        except:
            print ('ERROR: Negative value encountered at', time)
            sys.exit(1)"""
    print(compositions)
            
    return compositions #will return the last composition value of all particles that were a part of the sim no matter if they last the whole time

def write_output(compositions):
    f = open("/Users/nferich/GitHub/REBOUND_fragmentation/test_case/composition_output.txt", "w")
    for body in compositions:
        proper = '[%s]' % ' '.join(map(str, body))
        f.write(proper[1:-1]+"\n")
    f.close()

write_output(track_composition())


print("--- %s seconds ---" % (time.time() - start_time))
