import csv
import matplotlib.pyplot as plt
import numpy as np

Me_to_kg = 5.972e24
Re_to_m = 6.371e6

def error_propagation_division(x, y, z, delt_x, delt_y):
    
    propagated_error = z*((((delt_x/x)**2)+((delt_y/y)**2))**0.5)
    return(propagated_error)
    
    
def error_propagation_exponent(x, delt_x, exponent):
    
    propagated_error = delt_x*exponent*(x**(exponent-1))
    return(propagated_error)
    

data_file = "/Users/nferich/GitHub/REBOUND_fragmentation/exoplanet_archive_data/PSCompPars_2022.07.31_13.05.25.csv"

fields = []
raw_data = []

with open(data_file, 'r') as csvfile: 
    csv_reader = csv.reader(csvfile)
    
    line_count = 1
    
    for row in csv_reader:
        if line_count == 37:
            fields.append(row) #adds the list of data types for each column to the appropriate list
        elif line_count > 37:
            raw_data.append(row) #adds the data for each planet to the list
        line_count+=1

raw_data.pop(2299) #Kepler-1408 b - has a radius of .8900 R_e and an upper error of 68.91
           
deletion_list = []
for planet in raw_data:
    
    planet_mass = planet[6]
    planet_radius = planet[10]
    star_mass = planet[22]
    
    if int(planet[1]) > 1 or planet_mass == '' or planet_radius == '' or star_mass == '': #if the number of stars the planet orbits is greater than 1 or if it's missing a necessary measurement
        deletion_list.append(planet) #add that planet to the list of items to delete
    elif float(planet_radius) > 1.7: #IMPORTANT - Cutoff for when planet is considered terrestial or not
        deletion_list.append(planet) #add that planet to the list of items to delete
    elif float(planet_mass) > 2.0:
        deletion_list.append(planet) #add that planet to the list of items to delete
    else:
        continue


#removes incompatible data from list
for planet in deletion_list:
    raw_data.remove(planet)     
        
planet_radii = np.array([])
planet_radii_lower_unc = np.array([])
planet_radii_upper_unc = np.array([])       
planet_masses = np.array([])
planet_masses_upper_unc = np.array([])
planet_masses_lower_unc = np.array([])
star_masses = np.array([])
star_masses_upper_unc = np.array([])
star_masses_lower_unc = np.array([])

for planet in raw_data:
    
    planet_radii = np.append(planet_radii, float(planet[6]))
    if planet[7] == '': #if no error available
        planet_radii_upper_unc = np.append(planet_radii_upper_unc, 0) #just say it's 0
    else: #if there is an error available
        planet_radii_upper_unc = np.append(planet_radii_upper_unc, float(planet[7]))
    if planet[8] == '':
        planet_radii_lower_unc = np.append(planet_radii_lower_unc, 0)
    else:
        planet_radii_lower_unc = np.append(planet_radii_lower_unc, -1*float(planet[8]))
        
    
    planet_masses = np.append(planet_masses, float(planet[10]))
    if planet[11] == '':
        planet_masses_upper_unc = np.append(planet_masses_upper_unc, 0)
    else:
        planet_masses_upper_unc = np.append(planet_masses_upper_unc, float(planet[11]))
    if planet[12] == '':
        planet_masses_lower_unc = np.append(planet_masses_lower_unc, 0)
    else:
        planet_masses_lower_unc = np.append(planet_masses_lower_unc, -1*float(planet[12]))
        
    star_masses = np.append(star_masses, float(planet[22]))
    if planet[23] == '':
        star_masses_upper_unc = np.append(star_masses_upper_unc, 0)
    else:
        star_masses_upper_unc = np.append(star_masses_upper_unc, float(planet[23]))
    if planet[24] == '':
        star_masses_lower_unc = np.append(star_masses_lower_unc, 0)
    else:
        star_masses_lower_unc = np.append(star_masses_lower_unc, -1*float(planet[24]))
        
 
planet_densities = np.array([]) #it's assumed the planet is uniform and spherical
planet_densities_upper_unc = np.array([]) 
planet_densities_lower_unc = np.array([]) 
     

for i in range(planet_masses.size): 
    
    mass = planet_masses[i]*Me_to_kg #changes units of planetary mass to kg
    radius = planet_radii[i]*Re_to_m #changes units of planetary radius to m
    
    volume = (4/3)*np.pi*(radius**3)
    density = mass/volume
    planet_densities = np.append(planet_densities, density)
    
    mass_lower_unc = planet_masses_lower_unc[i]*Me_to_kg
    mass_upper_unc = planet_masses_upper_unc[i]*Me_to_kg
    radius_lower_unc = planet_radii_lower_unc[i]*Re_to_m
    radius_upper_unc = planet_radii_upper_unc[i]*Re_to_m
    
    volume_lower_unc = error_propagation_exponent(radius, radius_lower_unc, 3)
    volume_upper_unc = error_propagation_exponent(radius, radius_upper_unc, 3)
    
    density_lower_unc = error_propagation_division(mass, volume, density, mass_lower_unc, volume_lower_unc)
    density_upper_unc = error_propagation_division(mass, volume, density, mass_upper_unc, volume_upper_unc)
    
    planet_densities_lower_unc = np.append(planet_densities_lower_unc, density_lower_unc) 
    planet_densities_upper_unc = np.append(planet_densities_upper_unc, density_upper_unc)
    
    
planet_densities_unc = np.vstack((planet_densities_lower_unc, planet_densities_upper_unc)) #makes a single array filled with the lower and upper erroes for densities
star_masses_unc = np.vstack((star_masses_lower_unc, star_masses_upper_unc)) #makes a single array filled with the lower and upper erroes for stellar masses

#plt.figure(figsize=(6,6))
#plt.scatter(planet_densities, star_masses, 10)
#plt.errorbar(planet_densities, star_masses, xerr=planet_densities_unc,  fmt='o', markersize=3, lw=.5)
plt.errorbar(planet_densities, star_masses, xerr=planet_densities_unc, yerr=star_masses_unc, fmt='o', color='black', ecolor='tab:gray',markersize=3, elinewidth=.4, capsize=.8 )
#plt.errorbar(planet_radii, star_masses, xerr=planet_radii_upper_unc, fmt='o')
plt.xlabel('Planetary Bulk Density ($kg/m^{3}$)')
plt.ylabel('Stellar Mass ($M_{\odot}$)')

plt.savefig('density_vs_Mstar_errors.pdf', bbox_inches='tight', pad_inches=1.25)
#plt.savefig('density_vs_Mstar_errors.png', bbox_inches='tight', pad_inches=1.25, dpi=300)	

plt.show()  

"""for i in planet_densities_upper_unc:
    if i >10000:
        print(np.where(planet_densities_upper_unc == i))
        
print(raw_data[80])
print(planet_densities[80])
print(planet_densities_upper_unc[80])

a = .770*Me_to_kg
b= 1.43*Re_to_m
v = (4/3)*np.pi*b**3
print(a/v)"""
    
    
#print(next(csv_reader))
    


    
    
    
