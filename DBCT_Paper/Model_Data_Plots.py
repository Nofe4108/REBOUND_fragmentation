"""
Produces plots 5 through 8 in Ferich et al. (IN PREP). The code pulls 
data from the DBCT while it's running
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy import units as u
import time
import sys

start_time = time.time() #Timer to see how long running this code takes

#Lists needed to make the plots - dc stands for disruptive collision
dc_v_over_vesc = []
dc_target_masses = []
dc_proj_masses = []
dc_lr_masses = []
dc_lr_cmfs = []
dc_target_cmfs = []
dc_fragment_masses = []
dc_B_over_Rtarg = []
dc_ideal_ejecta_CMF = []
dc_actual_ejecta_CMF = []
dc_collision_types = [] # Blue is accretive, orange is erosive

# Constants
G = 39.478 #Gravitational Constant in AU/YR/Msun
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth

# Units
kg_per_m_cubed = u.kg/u.m**3
Msun_per_au_cubed = u.Msun/u.au**3

# Densities
core_density = ((7874.0*kg_per_m_cubed).to(Msun_per_au_cubed)) # Density of iron
mantle_density = ((3000.0*kg_per_m_cubed).to(Msun_per_au_cubed)) # Density of mantle material

def organize_compositions(init_compositions):
    """Data organizing function
    
    Takes in raw string data extracted from mantle stripping input file
    Creates a list filled with properties of each particle
    The properties include hash, mass, and CMF
    
    Parameters:
    init_compositions (list) -- raw particle data from input file

    Returns:
    compositions (list) -- nested list with properly formatted particle data
    """
    compositions = []
    for particle in init_compositions: 
        particle_data = [] 
        for i in range(len(particle)): 
            if i == 0: 
                particle_data.append(int(particle[i])) # The particle's REBOUND hash
            elif i < 3: 
                particle_data.append(float(particle[i])) # Adds its mass than CMF
        compositions.append(particle_data) 
    return(compositions)

def calc_core_radius(mass, core_frac, core_density):
    """Returns the radius of a differentiated object's spherical, uniform core"""
    core_mass = mass*core_frac
    core_radius = ((3*core_mass)/(4*np.pi*core_density))**(1/3)
    return(core_radius)

def calc_radius(mass, core_frac, mantle_density, core_radius):
    """Returns the outer radius of a differentiated object with a known core radius"""
    mantle_frac = 1-core_frac
    mantle_mass = mass*mantle_frac
    radius = (core_radius**3+((3*mantle_mass)/(4*np.pi*mantle_density)))**(1/3)
    return(radius)

def calc_ejecta_core_frac(collision_type, impact_parameter, min_ejecta_cmf, max_ejecta_cmf, targ_radius, targ_core_radius, proj_radius, proj_core_radius):
    """Mantle Stripping Model
    
    Calculates CMF in collisional ejecta based on impact parameters
    See Ferich et al. (in prep) for details of model
    
    Parameters:
    collision_type (int) -- can be either accretive (2) or erosive (3 & 4)
    impact_parameter (float) -- B in model equations from Ferich et al. (in prep)
    min_ejecta_cmf (float) -- minimum possible CMF of impact ejecta
    max_ejecta_cmf (float) -- maximum possible CMF of impact ejecta
    targ_radius (float) -- outer radius of target
    targ_core_radius (float) -- radius of target's core
    proj_radius (float) -- outer radius of projectile
    proj_core_radius (float) -- radius of projectile's core
    
    Returns:
    ejecta_cmf (float) -- The ideal fraction of core material in the collision ejecta
    """
    ejecta_cmf = 0.0
    if collision_type == 2:
        if proj_core_radius != 0:
            slope = (max_ejecta_cmf-min_ejecta_cmf)/(-2*proj_core_radius)
            max_impact_parameter = targ_radius+proj_core_radius 
            min_impact_parameter = targ_radius-proj_core_radius 
            if impact_parameter >= max_impact_parameter: # Oblique collisions
                ejecta_cmf = min_ejecta_cmf 
            elif max_impact_parameter > impact_parameter > min_impact_parameter: 
                ejecta_cmf = min_ejecta_cmf+(slope*(impact_parameter-max_impact_parameter)) 
            elif impact_parameter <= min_impact_parameter: # Head-on collisions 
                ejecta_cmf = max_ejecta_cmf 
    elif collision_type == 3 or collision_type == 4:
        if targ_core_radius != 0:
            slope = (max_ejecta_cmf-min_ejecta_cmf)/(-2*targ_core_radius)
            max_impact_parameter = proj_radius+targ_core_radius
            min_impact_parameter = proj_radius-targ_core_radius 
            if impact_parameter >= max_impact_parameter:  # Oblique collisions
                ejecta_cmf += min_ejecta_cmf 
            elif max_impact_parameter > impact_parameter > min_impact_parameter: 
                ejecta_cmf += min_ejecta_cmf+(slope*(impact_parameter-max_impact_parameter)) 
            elif impact_parameter <= min_impact_parameter: # Head-on collisions
                ejecta_cmf += max_ejecta_cmf         
    return(ejecta_cmf)

def calc_ejecta_core_mantle_masses(targ_mass, lr_mass, obj_mass, obj_radius, obj_cmf, ejecta_cmf):
    """Calculates amount of core and mantle material in ejecta
    
    If the core or mantle of the stripped object doesn't have enough
    material to reach the ideal CMF for the collision ejecta, then
    this function will pull more material out of the other layer.
    Basically makes sure that the code doesn't put core or mantle material that
    the stripped object doesn't have into the ejecta.
    
    
    Parameters:
    targ_mass (float) -- Mass of the target
    lr_mass (float) -- Mass of the largest remnant
    obj_mass (float) -- Radius of the stripped object in collision
    obj_radius (float) -- Radius of the stripped object in collision
    obj_cmf (float) -- Core mass fraction of stripped object
    ejecta_cmf (float) -- the ideal CMF of the collision ejecta

    Returns:
    actual_ejecta_masses (list) -- non-ideal core and mantle masses of the ejecta
    """
    ejecta_mass = abs(targ_mass-lr_mass)
    ideal_ejecta_masses = [(ejecta_mass*ejecta_cmf),(ejecta_mass*(1-ejecta_cmf))] # Ideal masses of the ejecta's core and mantle
    obj_layer_masses = [obj_mass*obj_cmf, obj_mass*(1-obj_cmf)] # How much mass the object actually has in each layer
    actual_ejecta_masses = [] # Actual masses of the ejecta's core and mantle
    if obj_layer_masses[0] >= ideal_ejecta_masses[0]:
        actual_ejecta_masses.append(ideal_ejecta_masses[0])
    else:
        actual_ejecta_masses.append(obj_layer_masses[0])
        ideal_ejecta_masses[1] += (ideal_ejecta_masses[0]-obj_layer_masses[0])
    if obj_layer_masses[1] >= ideal_ejecta_masses[1]:
        actual_ejecta_masses.append(ideal_ejecta_masses[1])
    else:
        actual_ejecta_masses.append(obj_layer_masses[1])
        actual_ejecta_masses[0] += (ideal_ejecta_masses[1]-obj_layer_masses[1])
    return(actual_ejecta_masses)



def track_composition(collision_report_file, composition_input_file, ejection_file, impact_parameter_file, velocity_file, min_core_collision_frac, max_core_collision_frac):
    """Main function
    
    Tracks how collisions change the CMFs of objects from REBOUND fragmentation sim
    Returns the final compositions of all remaining objects 
    Has 2 main sections for mergers and disruptive collisions (accretive and erosive)
    
    
    Parameters:
    collision_report_file (str) -- pathway to collision report file
    collision_input_file (str) -- pathway to the mantle stripping input file
    ejection_file (str) -- pathway to file that lists objects ejected from sim
    impact_parameter_file (str) -- pathway to file that lists impact parameters of collisions
    velocity_file (str) -- pathway to file that lists velocities of each collision
    min_core_collision_frac (float) -- minimum possible CMF of impact ejecta
    max_core_collision_frac (float) -- maximum possible CMF of impact ejecta

    Returns:
    compositions (list) -- nested list with compositional data of final objects
    """
    # Extracts compositional data
    f = open(composition_input_file, 'r')
    raw_compositions = [line.split() for line in f.readlines()]
    compositions = organize_compositions(raw_compositions) # This list keeps track of all the particle data
    for obj in compositions:
        if obj[2] > 1.0 or obj[2] < 0.0:
            print ('ERROR: CMF does not have a realistic value')
            sys.exit(1) 
    f.close() 
    
    # Extracts collisional data  
    f = open(collision_report_file, 'r')
    collision_blocks = f.read().split("\n")
    collisions = [block for block in collision_blocks if len(block) > 0]
    f.close()
    
    # Extracts the impact parameters for each collision
    f = open(impact_parameter_file, 'r')
    b = f.read() #reads each line in file and splits its value into an array - these impact parameters are divided by R_targ
    impact_parameters_raw = b.split("\n")
    impact_parameters = [x.strip('b/Rtarg:     ') for x in impact_parameters_raw]
    f.close()
    
    # Extracts the velocties of each collision - only used for plotting
    f = open(velocity_file, 'r')
    v = f.read() #reads each line in file and splits its value into an array - these impact parameters are divided by R_targ
    velocities_raw = v.split("\n")
    collision_velocities = [x.strip('Vimp/Vesc:     ') for x in velocities_raw]
    
    destroyed_object_hashes = [] # Keeps track of objects that get destroyed in collisions
    
    ############### START OF MAIN LOOP ############### 
    # Iterates through every collision from sim
    # Determines the change in objects' compositions from each collision
    # Adds new fragments to the compositions list when necessary
    for i in range(len(collisions)): 
        collision = collisions[i].split() 
        time = float(collision[0]) 
        collision_type = int(collision[1]) 
        if collision_type == 0: # Elastic bounces do nothing so loop moves to next collision
            continue
        target_hash = int(collision[2]) 
        largest_remnant_mass = float(collision[3]) # Mass of the target after collision
        proj_hash = int(collision[4]) 
        no_frags = int((len(collision)-5)/2) #gives number of frags formed in collision - gets rid of 5 values from length which are basically about the collision type and target and projectile then divides in half since the list contains both the fragment hash and its mass
        frag_hashes = [int(collision[j*2+3]) for j in range(1,no_frags+1)] #list of the hashes of the fragments - jumps from hash to hash using (i*2+3) - range length uses no_frags variable to get the correct amount of hashes
        frag_masses = [float(collision[j*2+4]) for j in range(1,no_frags+1)]

        # Determines if object collided with something not in the compositions list (Star, Jupiter, etc.)
        # If so, just destroys projectile and moves to the next collision
        # If fragments created, adds those to compositions list with CMF of 0.3
        #FIXME: Check to make sure this doesn't break anything
        big_obj_collision_flag = 1
        for j in range(len(compositions)):
            if target_hash == compositions[j][0]:
                big_obj_collision_flag += -1
                break
        if big_obj_collision_flag == 1:
            destroyed_object_hashes.append(proj_hash)
            if no_frags != 0:
                for k in range(no_frags):
                     frag_data = [frag_hashes[k], frag_masses[k], 0.3]
                     compositions.append(frag_data) 
            continue   

        targ_idx = [idx for idx in range(len(compositions)) if int(compositions[idx][0])==target_hash][0] # Index of target in the compositions list
        proj_idx = [idx for idx in range(len(compositions)) if int(compositions[idx][0])==proj_hash][0] # Index of projectile in the compositions list
        target_mass = float(compositions[targ_idx][1]) 
        proj_mass = float(compositions[proj_idx][1]) 
        target_core_frac = compositions[targ_idx][2]
        proj_core_frac = compositions[proj_idx][2] 
        
        # This sequence estimates the core radii of the target and projectile and the impact parameter
        target_sim_radius = calc_core_radius(target_mass, 1.0, mantle_density)
        proj_sim_radius = calc_core_radius(proj_mass, 1.0, mantle_density)
        target_core_radius = calc_core_radius(target_mass, target_core_frac, core_density)
        target_radius = calc_radius(target_mass, target_core_frac, mantle_density, target_core_radius)
        proj_core_radius = calc_core_radius(proj_mass, proj_core_frac, core_density)
        proj_radius = calc_radius(proj_mass, proj_core_frac, mantle_density, proj_core_radius)
        sim_impact_param = float(impact_parameters[i])*target_sim_radius
        sine_of_impact_angle = sim_impact_param/(target_sim_radius+proj_sim_radius)
        impact_parameter = (target_radius+proj_radius)*sine_of_impact_angle
                                                
    ############### PERFECT MERGER ############### 
        if collision_type == 1: #perfect merger
            compositions[targ_idx][2] = ((target_core_frac*target_mass)+(proj_core_frac*proj_mass))/largest_remnant_mass #changes the composition fraction for each specie in the target - basically weighted average of initial target compoisition and mass with the projectile composition and mass
 
    ############### DISRUPTIVE COLLISIONS ############### 
        elif collision_type == 2 or collision_type == 3 or collision_type == 4:
            ejecta_core_frac = calc_ejecta_core_frac(collision_type, impact_parameter, min_core_collision_frac, max_core_collision_frac, target_radius, target_core_radius, proj_radius, proj_core_radius)
            
            # Collects impact data for Disruptive Collisions
            if target_mass != largest_remnant_mass and sum(frag_masses) != proj_mass:
                dc_v_over_vesc.append(float(collision_velocities[i]))
                dc_target_masses.append(target_mass*mass_conversion)
                dc_proj_masses.append(proj_mass*mass_conversion)
                dc_lr_masses.append(largest_remnant_mass*mass_conversion)
                dc_target_cmfs.append(target_core_frac)
                dc_fragment_masses.append([frag_masses[j]*mass_conversion for j in range(len(frag_masses))])
                dc_B_over_Rtarg.append(impact_parameter/target_radius)
                dc_ideal_ejecta_CMF.append(ejecta_core_frac)
            
            if collision_type == 2:
                ejecta_layer_masses = calc_ejecta_core_mantle_masses(target_mass, largest_remnant_mass, proj_mass, proj_radius, proj_core_frac, ejecta_core_frac)
                compositions[targ_idx][2] = ((target_mass*target_core_frac)+(ejecta_layer_masses[0]))/largest_remnant_mass # Changes target CMF to largest remnant CMF for accretion
                
                # Collects impact data for Accretive Collisions
                if target_mass != largest_remnant_mass and sum(frag_masses) != proj_mass:
                    dc_collision_types.append("tab:blue")
                    dc_lr_cmfs.append(compositions[targ_idx][2])
                    if abs(target_mass-largest_remnant_mass) == 0:
                        dc_actual_ejecta_CMF.append(0)
                    else:
                        dc_actual_ejecta_CMF.append(ejecta_layer_masses[0]/abs(target_mass-largest_remnant_mass))
                
            else:
                ejecta_layer_masses = calc_ejecta_core_mantle_masses(target_mass, largest_remnant_mass, target_mass, target_radius, target_core_frac, ejecta_core_frac)
                compositions[targ_idx][2] = ((target_mass*target_core_frac)-(ejecta_layer_masses[0]))/largest_remnant_mass # Changes target CMF to largest remnant CMF for erosion
                
                # Collects impact data for Erosive Collisions
                if target_mass != largest_remnant_mass and sum(frag_masses) != proj_mass:
                    dc_collision_types.append("tab:orange")
                    dc_lr_cmfs.append(compositions[targ_idx][2])
                    if abs(target_mass-largest_remnant_mass) == 0:
                        dc_actual_ejecta_CMF.append(0)
                    else:
                        dc_actual_ejecta_CMF.append(ejecta_layer_masses[0]/abs(target_mass-largest_remnant_mass))
            
            # Target CMF calculation tolerance
            if compositions[targ_idx][2] - 1.0 > 0.0 and compositions[targ_idx][2] - 1.0 < 1.0e-5: # If target CMF is just above 1.0
                compositions[targ_idx][2] = 1.0
            if compositions[targ_idx][2] < 0.0 and compositions[targ_idx][2] > -1.0e-5: # If target CMF is just below 0.0
                compositions[targ_idx][2] = 0.0
            
            total_core_mass = (target_core_frac*target_mass)+(proj_core_frac*proj_mass) 
            largest_remnant_core_mass = largest_remnant_mass*compositions[targ_idx][2] 
            total_frag_core_mass = total_core_mass - largest_remnant_core_mass 
            total_frag_mass = (target_mass+proj_mass) - largest_remnant_mass 
            frag_core_frac = total_frag_core_mass/total_frag_mass # All fragments get the same CMF
            
            # Fragment CMF calculation tolerance
            if frag_core_frac - 1.0 > 0.0 and frag_core_frac - 1.0 < 1.0e-5: # If fragment CMF is just above 1.0
                frag_core_frac += 1.0 - frag_core_frac
            if frag_core_frac < 0.0 and frag_core_frac > -1.0e-5: # If fragment CMF is just below 0.0
                frag_core_frac += 0.0 - frag_core_frac
            
            # Fragment CMF error
            if frag_core_frac < 0.0:
                print ('ERROR: Fragment CMF is negative at time: ', time)
                sys.exit(1)
            elif frag_core_frac > 1.0:
                print ('ERROR: Fragment CMF is greater than 1.0 at time: ', time)
                sys.exit(1)
                
            for j in range(no_frags):
                if frag_masses[j] < 0:
                    print ('ERROR: Fragment mass is negative at time: ', time)
                    sys.exit(1)
                else:
                    frag_data = [frag_hashes[j], frag_masses[j], frag_core_frac]
                    compositions.append(frag_data) # Fragment now added to the object tracking list
                
        compositions[targ_idx][1] = largest_remnant_mass # Changes target mass to largest remnant mass
        
        # Checks to see if there are any errors with the largest remnant's properties 
        for j in range(len(compositions[targ_idx])):
            if j == 1:
                if compositions[targ_idx][j] < 0.0: 
                    print ('ERROR: Largest remnant mass is negative at time: ', time)
                    sys.exit(1)
            elif j > 1:
                if compositions[targ_idx][j] < 0.0:
                    print ('ERROR: Largest remnant CMF is negative at time: ', time)
                    sys.exit(1)
                elif compositions[targ_idx][j] > 1.0:
                    print ('ERROR: Largest remnant CMF is greater than 1.0 at time: ', time)
                    sys.exit(1)
                    
                    
        # Checks to make sure the projectile isn't a second largest remnant before deletion
        for hsh in frag_hashes:
            if hsh != proj_hash and hsh == frag_hashes[-1]:
                destroyed_object_hashes.append(proj_hash) 
            elif hsh == proj_hash:
                break
            else:
                continue
        
        if collision_type == 1:  
            destroyed_object_hashes.append(proj_hash) # The projectile in a merger is always destroyed
    
    ############### END OF MAIN LOOP ###############
    
    # Removes destroyed objects from compositions list    
    for hsh in destroyed_object_hashes:
        for i in range(len(compositions)):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break
            else:
                continue
    
    # Removes ejected objects from compositions list
    f = open(ejection_file, 'r')
    ejections_raw = [line.split() for line in f.readlines()]
    ejections = [abs(int(ejections_raw[j][2])) for j in range(len(ejections_raw))]
    for hsh in ejections:
        for i in range(len(compositions)):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break
            else:
                continue
    f.close()
         
    return(compositions)

######### DATA COLLECTION #################

min_cmf = 0.0 #mimimum fraction of core material in ejecta
max_cmf = 1.0 #maximum fraction of core material in ejecta

# File pathways
collision_file_pw = "new_collision_reports/new_collision_report"
composition_input_file_pw = "DBCT_input/uni_DBCT_input"
impacter_param_file_pw = "impact_parameters/impact_parameters"
ejec_file_pw = "ejections/ejections"
velocity_file_pw = "velocities/velocities"

file_range = np.arange(1,51,1) # Number of simulation files used for the plots

for i in file_range:
    collision_file = collision_file_pw + str(i) + ".txt"
    composition_input_file = composition_input_file_pw + ".txt"
    impacter_param_file = impacter_param_file_pw + str(i) + ".txt"
    ejec_file = ejec_file_pw + str(i) + ".txt"
    velo_file = velocity_file_pw + str(i) + ".txt"
    final_composition = track_composition(collision_file, composition_input_file, ejec_file, impacter_param_file, velo_file, min_cmf, max_cmf) 

# This sequence is used to calculate the impact velocity in km/sec
total_radius = [calc_radius(dc_target_masses[i]/mass_conversion, 0, mantle_density, 0)+calc_radius(dc_proj_masses[i]/mass_conversion, 0, mantle_density, 0) for i in range(len(dc_target_masses))] #supposed to have units of au
total_mass = [(dc_target_masses[i]+dc_proj_masses[i])/mass_conversion for i in range(len(dc_target_masses))]    
v_esc_au_yr = [np.sqrt(2*G*total_mass[i]/total_radius[i]) for i in range(len(dc_target_masses))]
v_esc = [v_esc_au_yr[i]*(u.au/u.yr).to(u.km/u.s) for i in range(len(dc_target_masses))] # Should be km/s
v_impact = [dc_v_over_vesc[i]*v_esc[i] for i in range(len(dc_target_masses))]

######### PLOTTING ##################

plt.rcParams['axes.axisbelow'] = True # Makes sure grid is behind points 

# Figure 5 
fig1, (ax1_1, ax1_2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6,10))
fig1.subplots_adjust(hspace=0.03)
ax1_1.scatter(dc_B_over_Rtarg, dc_ideal_ejecta_CMF, color=dc_collision_types, s=5.0, marker='o', linewidths=0, alpha=0.6)
ax1_1.set_ylabel('Ideal Ejecta CMF', fontsize='large') 
ax1_1.grid(alpha=0.7)  
ax1_1.tick_params(length=3, width=1, which='major')
ax1_1.minorticks_on()
ax1_1.set_xlim([-0.05, 2.05])
ax1_1_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='Accretive Collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='Erosive Collision', markerfacecolor='tab:orange', markersize=10.0)]
legend1 = ax1_1.legend(handles=ax1_1_legend_elements, loc = 'upper right', framealpha = .7)
ax1_2.scatter(dc_B_over_Rtarg, dc_actual_ejecta_CMF, color=dc_collision_types, s=5.0, marker='o', linewidths=0, alpha=0.6)
ax1_2.set_xlabel('B/$R_{targ}$', fontsize='large')
ax1_2.set_ylabel('Actual Ejecta CMF', fontsize='large')  
ax1_2.grid(alpha=0.7) 
ax1_2.minorticks_on()
ax1_2.set_xlim([-0.05, 2.05])
plt.savefig("graphs/ejecta_cmfs_vs_impact_parameter.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/ejecta_cmfs_vs_impact_parameter.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/ejecta_cmfs_vs_impact_parameter.png', bbox_inches='tight', dpi=300)

# Figure 6
fig2, ax2 = plt.subplots(figsize=(6,5))
ax2.scatter(dc_target_masses, dc_v_over_vesc, c=dc_collision_types, s=5.0, marker='o', linewidths=0, alpha=0.6, zorder=2)
ax2.set_xlabel('Target Mass ($M_{\u2295}$)', fontsize='large')
ax2.set_ylabel('$v_{imp}/v_{esc}$', fontsize='large')
ax2.grid(alpha=0.7)
plt.xscale("log")  
plt.yscale("log")
ax2.axvline(.093, label='Embryo Mass', color='black', linestyle=':', alpha=0.5, zorder=1)
ax2.axvline(.0093, label = 'Planetesimal Mass', color = 'black', linestyle = '-.', alpha=0.5, zorder=1)
ax2.axvline(0.00465, label = 'Minimum Fragment Mass', color = 'black', linestyle = '--', alpha=0.5, zorder=1)
ax2.set_ylim([0.95,105.1])
ax2_legend_elements_1 = [Line2D([], [], color='black', linestyle=':',label='Embryo Mass', alpha=0.5),
                      Line2D([], [], color='black', linestyle='-.', label='Planetesimal Mass', alpha=0.5),
                      Line2D([], [], color='black', linestyle='--', label='Minimum Fragment Mass', alpha=0.5),
                      Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='Accretive Collision', markerfacecolor='tab:blue', markersize=7.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='Erosive Collision', markerfacecolor='tab:orange', markersize=7.0)] 
"""ax2_legend_elements_2 = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='Accretive Collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='Erosive Collision', markerfacecolor='tab:orange', markersize=10.0)] """
legend2_1 = plt.legend(handles=ax2_legend_elements_1, loc='upper right', framealpha=0.7, fontsize=8.0)
#Legend2_2 = plt.legend(handles=ax2_legend_elements_2, loc = 'upper left', framealpha = 0.7, fontsize=8.0)
plt.gca().add_artist(legend2_1)
#plt.gca().add_artist(Legend2_2)
plt.savefig("graphs/target_mass_vs_normalied_vimp.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/target_mass_vs_normalied_vimp.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/target_mass_vs_normalied_vimp.png', bbox_inches='tight', dpi=300)

"""# Figure 7
fig3, (ax3_1, ax3_2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6,10))
fig3.subplots_adjust(hspace=0.03)
ax3_1.scatter(dc_v_over_vesc, dc_ideal_ejecta_CMF, c=dc_collision_types, s=5.0, marker='o', linewidths=0, alpha=0.7)
ax3_1.set_ylabel('Ideal Ejecta CMF', fontsize='large')
plt.xscale("log")
ax3_1.minorticks_on()
ax3_1.set_xlim([0.95, 105.1])
ax3_1.grid(alpha=0.7)
ax3_2.scatter(dc_v_over_vesc, dc_actual_ejecta_CMF, c=dc_collision_types, s=5.0, marker='o', linewidths=0, alpha=0.7)
ax3_2.set_xlabel('$v_{imp}/v_{esc}$', fontsize='large')
ax3_2.set_ylabel('Actual Ejecta CMF', fontsize='large')
plt.xscale("log")
ax3_2.set_xlim([0.95, 105.1])
ax3_2.grid(alpha=0.7)
ax3_2_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='Accretive collision', markerfacecolor='tab:blue', markersize=5.5),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='Erosive collision', markerfacecolor='tab:orange', markersize=5.5)] 
legend3 = ax3_2.legend(handles=ax3_2_legend_elements, loc = 'upper left', framealpha = 0.4, fontsize=6.5)
plt.savefig("graphs/fig7.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig7.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig7.png', bbox_inches='tight', dpi=300)

# Figure 8
fig4, ax4 = plt.subplots(figsize=(6,5))
ax4.scatter(dc_target_cmfs, dc_lr_cmfs, c=dc_collision_types, s=3.0, marker='o', linewidths=0, alpha=0.6)
ax4.set_xlabel('Target CMF', fontsize='large')
ax4.set_ylabel('Largest Remnant CMF', fontsize='large')
ax4.set_xlim([-0.01, 1.01])
ax4.set_ylim([-0.01, 1.01])
ax4.grid(alpha=0.7)
ax4.minorticks_on()
ax4_legend_elements = [Line2D([], [], color='tab:blue', marker='o', lw=0.0, label='accretive collision', markerfacecolor='tab:blue', markersize=10.0),
                      Line2D([], [], color='tab:orange', marker='o', lw=0.0, label='erosive collision', markerfacecolor='tab:orange', markersize=10.0)] 
legend = plt.legend(handles=ax4_legend_elements, loc = 'upper right', framealpha = .7)
plt.savefig("graphs/fig8.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig8.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/fig8.png', bbox_inches='tight', dpi=300)

print("--- %s seconds ---" % (time.time() - start_time))"""