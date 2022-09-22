import numpy as np
from numpy import random


f = open("/Users/nferich/GitHub/REBOUND_fragmentation/disk_simulations/disk_mantle_stripping_input.txt", "w")

no_bodies = 5 #number of bodies in disk simulation
mass = 1.0

for i in range(no_bodies): #goes through each particle in compositions

    core_frac = random.randint(0, 0.4e3)/1e3
    mantle_frac = 1-core_frac
    
    f.write(str(i) + ' ')
    f.write(str(mass) + ' ')
    f.write(str(core_frac) + ' ')
    f.write(str(mantle_frac) + ' ')
    f.write('\n')#go to a new line and to a new particle
f.close()