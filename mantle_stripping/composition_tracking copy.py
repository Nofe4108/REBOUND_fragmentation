#!/usr/bin/env python

import numpy as np
import time
import sys

start_time = time.time()

def organize_compositions(init_compositions):
    compositions = []
    for line in init_compositions: #iterates through each line from file
        particle_data = [] #This will hold data for individual particle - will contain its hash, mass, and the fractions for each layer
        for i in range(len(line)): #goes through the index for each element in a line
            if i == 0: #First value in each line is the hash so that'll be an integer
                particle_data.append(int(line[i]))
            elif i < 3: #The other values will be floats - doesn't add initial semi-major axis or eccentricity to list
                particle_data.append(float(line[i]))
        compositions.append(particle_data) #add all the data about the particle to the compositions list
    return(compositions)

def track_composition(comp_input_file, collision_file, ejection_file):
    destroyed_object_hashes = []
    f = open(comp_input_file, 'r')
    init_compositions = [line.split() for line in f.readlines()]
    """for line in init_compositions:
        try:
            sum([float(x) for x in line[2:]]) == 1.0
        except:
            print ('ERROR: Realative abundances do not add up to 1.0')
            sys.exit(1)
    try:
        init_hashes = [int(x[0]) for x in init_compositions]
    except:
        init_hashes = [x[0].value for x in init_compositions]"""
    no_species = len(init_compositions[0])-4
    file = open(collision_file, 'r')
    blocks = file.read().split("\n")
    blocks = [block for block in blocks if len(block) > 0]
    compositions = organize_compositions(init_compositions)
    print(compositions[0])
    for i in range(len(blocks)):
        block = blocks[i].split()
        time = float(block[0])
        collision_type = int(block[1])
        if collision_type == 0:
            continue
        target_hash = int(block[2])
        target_mass = float(block[3])
        projectile_hash = int(block[4])
        
        big_object_collision = 1
        for i in range(len(compositions)):
            if target_hash==compositions[i][0]:
                big_object_collision+= -1
                break
        if big_object_collision==1:
            proj_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==projectile_hash][0]  #same thing as above except for projectile
            proj_mass = float(compositions[proj_idx][1]) #same thing as above except for projectile
            proj_core_frac = compositions[proj_idx][2] #same thing as above except for projectile
            destroyed_object_hashes.append(projectile_hash)
            continue 

        targ_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==target_hash][0]
        proj_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==projectile_hash][0]    
        last_target_mass = float(compositions[targ_idx][1])
        last_proj_mass = float(compositions[proj_idx][1])
        last_target_abundances = compositions[targ_idx][-no_species:]   
        last_projectile_abundances = compositions[proj_idx][-no_species:] 
        no_frags = int((len(block)-5)/2)
        frag_hashes = [int(block[i*2+3]) for i in range(1,no_frags+1)]
        frag_masses = [float(block[i*2+4]) for i in range(1,no_frags+1)]
        if collision_type == 1: #perfect merger
            for i in range(no_species):
                compositions[targ_idx][i+2]=(float(last_target_abundances[i])*last_target_mass+float(last_projectile_abundances[i])*last_proj_mass)/target_mass
        if collision_type == 2: #partial accretion
            mass_accreted = target_mass-last_target_mass
            for i in range(no_species):
                compositions[targ_idx][i+2]=(float(last_target_abundances[i])*last_target_mass+float(last_projectile_abundances[i])*mass_accreted)/target_mass
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j]]+last_projectile_abundances
                compositions.append(frag_data)
                try:
                     any(n < 0 for n in frag_data) == False
                except:
                    print ('ERROR: Negative value encountered in frag data at', time)
                    sys.exit(1)
        if collision_type == 3 or collision_type == 4: #partial errosion, target abundances stay the same
            mass_lost = last_target_mass-target_mass
            frag_abundances = [(float(last_target_abundances[i])*mass_lost+float(last_projectile_abundances[i])*last_proj_mass)/np.sum(frag_masses) for i in range(no_species)]
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j]]+frag_abundances
                try:
                     any(n < 0 for n in frag_data) == False
                except:
                    print ('ERROR: Negative value encountered in frag data at', time)
                    sys.exit(1)
                compositions.append(frag_data)
        compositions[targ_idx][1]=target_mass
        try: 
            any(n < 0 for n in compositions[targ_idx]) == False
        except:
            print ('ERROR: Negative value encountered at', time)
            sys.exit(1)
            
        #checks to make sure projectile isn't a second largest remnant before deleting it
        for hsh in frag_hashes:
            if hsh != projectile_hash and hsh == frag_hashes[-1]:
                destroyed_object_hashes.append(projectile_hash) #IMPORTANT CHANGE - removes projectile particle from the compositions list if it gets destroyed in collision  
            elif hsh == projectile_hash:
                break
            else:
                continue
        if collision_type == 1: #Does the same thing as loop above - specifically for a merger 
            destroyed_object_hashes.append(projectile_hash) #IMPORTANT CHANGE - removes projectile particle from the compositions list if it gets destroyed in collision


    #Destroys particles by hash      
    for hsh in destroyed_object_hashes:
        for i in range(len(compositions)):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break
            else:
                continue

    #ADDED THIS TO REMOVE EJCTED OBJECTS
    e_file = open(ejection_file, 'r')
    ejections_raw = [line.split() for line in e_file.readlines()]
    ejections = [int(ejections_raw[i][2]) for i in range(len(ejections_raw))]
    
    for hsh in ejections:
        for i in range(len(compositions)):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break
            else:
                continue
  
    e_file.close()
    print(len(compositions))
    
    return compositions

def write_output( compositions, output_file):
    f = open(output_file, "w")
    for body in compositions:
        proper = '[%s]' % ' '.join(map(str, body))
        f.write(proper[1:-1]+"\n")
    f.close()

collision_file_pw = "/Users/nferich/Desktop/anna_collision_reports/new_collision_reports/new_collision_report"
comp_input_file_pw = "/Users/nferich/Desktop/anna_collision_reports/mantle_stripping_input/exp_mantle_stripping_input"
ejec_file_pw = "ejections/ejections"
output_file_pw = "comp_tracking_output/exp_comp_tracking_output"



file_range = np.arange(1,51,1)

for i in file_range:
    collision_file = collision_file_pw + str(i) + ".txt"
    comp_input_file = comp_input_file_pw + str(1) + ".txt"
    ejec_file = ejec_file_pw + str(i) + ".txt"
    output_file = output_file_pw + str(i) + ".txt"
    write_output(track_composition(comp_input_file, collision_file, ejec_file), output_file)


print("--- %s seconds ---" % (time.time() - start_time))
