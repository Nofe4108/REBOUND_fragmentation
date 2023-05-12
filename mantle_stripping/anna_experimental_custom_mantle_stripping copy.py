import numpy as np
import math
import time
import sys

#hello = [x.strip(' ') for x in hello]

#SIM CONVERSIONS
distance_conversion = 1.49597870691e11 #m to au
mass_conversion = 1.9885e30 #kg to Msun
time_conversion = 31557600 #sec to yrs

core_density = 7874.0*(distance_conversion**3/mass_conversion) #kg m^-3 - density of iron
mantle_density = 3000.0*(distance_conversion**3/mass_conversion) #kg m^3 -  value used in the simulation

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
def track_composition(min_frac, max_frac): #Main function that gives the compositions that will be outputted to the file
    f = open("/Users/nferich/Desktop/anna_collision_reports/inputs/disk_mantle_stripping_input1.txt", 'r')
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
    
     
    file = open("collision_report1.txt", 'r')
    blocks = file.read().split("\n")#pulls all the data out of the collision report - the list element for each collision is one big string 
    blocks = [block for block in blocks if len(block) > 0] #gets rid of the empty string at the end of the list
    
    #Gets the impact parameters for each collision
    b_file = open("output/impact_parameters1.txt", 'r')
    b = b_file.read() #reads each line in file and splits its value into an array - these impact parameters are divided by R_targ
    impact_parameters_raw = b.split("\n")
    impact_parameters = [x.strip('b/Rtarg:     ') for x in impact_parameters_raw]

    destroyed_object_hashes = [] #list of objects that get destroyed in a collision
    
    #START OF BIG LOOP
    for i in range(len(blocks)): #iterates through each value in blocks list - THIS IS A VERY BIG LOOP THAT CONTAINS ALL OF THE MASS TRANSFER DECISION MAKING
        block = blocks[i].split() #seperates each long string full of the collision data into its own list to be parsed through
        time = float(block[0]) #time of collision is first value
        collision_type = int(block[1]) #type of collision is second
        if collision_type == 0: #This is just an elastic bounce so it just continues to the next set of collision data in the blocks list
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
        
                
            compositions[targ_idx][2] = ((target_core_frac*target_mass)+(proj_core_frac*proj_mass))/largest_remnant_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile composition and mass
            compositions[targ_idx][3] = 1 - compositions[targ_idx][2]
            
 
 ####################### PARTIAL ACCRETION ######################
   
        if collision_type == 2: #partial accretion
        
            mass_accreted = largest_remnant_mass-target_mass #change in mass of the target after the collision - this time mass is added to target
 
            min_core_collision_frac = min_frac #minimum fraction of core material largest remnant gained in ejecta
            max_core_collision_frac = max_frac #maximum fraction of core material largest remnant gained in ejecta
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
            
            #If the target_core_frac is just barely above 1 or just barely below 0, it will be rounded down to 1 or up to 0
            if compositions[targ_idx][2] - 1 > 0 and compositions[targ_idx][2] - 1 < 1e-8: 
                compositions[targ_idx][2] = 1.0
                compositions[targ_idx][3] = 0.0
            if compositions[targ_idx][2] < 0 and compositions[targ_idx][2] > -1e-8:
                compositions[targ_idx][2] = 0.0
                compositions[targ_idx][3] = 1.0
            
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
            
            min_core_collision_frac = min_frac #minimum fraction of core material largest remnant lost in ejecta
            max_core_collision_frac = max_frac #maximum fraction of core material largest remnant lost in ejecta
            ejecta_core_frac = 0 #fraction of the eroded mass that is composed of core material (will depend on impact parameter)
            
            if target_core_radius == 0: #makes sure the core radius of the target isn't 0
                ejecta_core_frac = 0
            else:
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
            
            #If the target_core_frac is just barely above 1 or just barely below 0, it will be rounded down to 1 or up to 0
            if compositions[targ_idx][2] - 1 > 0 and compositions[targ_idx][2] - 1 < 1e-8: 
                compositions[targ_idx][2] = 1.0
                compositions[targ_idx][3] = 0.0
            if compositions[targ_idx][2] < 0 and compositions[targ_idx][2] > -1e-8:
                compositions[targ_idx][2] = 0.0
                compositions[targ_idx][3] = 1.0
            
            #If the frag_core_frac is just barely above 1 or just barely below 0, it will be rounded down to 1 or up to 0
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
        for i in range(len(compositions[targ_idx])):
            if i == 1:
                if compositions[targ_idx][i] < 0: 
                    print ('ERROR: Negative value for target mass encountered at', time)
                    sys.exit(1)
            elif i > 1:
                if compositions[targ_idx][i] < 0:
                    print ('ERROR: Negative value encountered in target data at', time)
                    sys.exit(1)
    
    #Destroys particles by hash      
    for hsh in destroyed_object_hashes:
        for i in range(len(compositions)):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break
            else:
                continue
    
    e_file = open("ejections/ejections1.txt", 'r')
    ejections_raw = [line.split() for line in e_file.readlines()]
    ejections = [int(ejections_raw[i][2]) for i in range(len(ejections_raw))]
    for hsh in ejections:
        for i in range(len(compositions)):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break
            else:
                continue
    
    
    return(compositions)

######## DATA COLLECTING LOOP ###########
        
"""final_compositions = []
min_core_collision_fracs = np.arange(0, 1.1, 0.1)
max_core_collision_fracs = np.arange(0, 1.1, 0.1)
for i in min_core_collision_fracs:
    fracs = []
    for j in max_core_collision_fracs:
        if j >= i: #ensures the max collision frac is never less than the minimum
            final_composition = track_composition(i, j)
            final_core_frac = sum([final_composition[k][2]/len(final_composition) for k in range(len(final_composition))])
            fracs.append(final_core_frac)
    final_compositions.append(fracs)"""
    
final_compositions = []
min_core_collision_fracs = np.arange(0, 1.1, 0.1)
max_core_collision_fracs = np.arange(0, 1.1, 0.1)
for i in min_core_collision_fracs:
    fracs = []
    for j in max_core_collision_fracs:
        if j >= i: #ensures the max collision frac is never less than the minimum
            final_composition = track_composition(i, j)
            final_embryo_composition = []
            for obj in final_composition:
                if obj[1]*330000 > 0.093: #if the object mass is bigger than original embryo mass
                    final_embryo_composition.append(obj)
            final_core_frac = sum([final_embryo_composition[k][2]/len(final_embryo_composition) for k in range(len(final_embryo_composition))])
            fracs.append(final_core_frac)
    final_compositions.append(fracs)

print(final_compositions) 

        

print("--- %s seconds ---" % (time.time() - start_time))
