#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 17:04:21 2022

@author: nferich
"""
f = open("/Users/nferich/Desktop/anna_collision_reports/REBOUND_output/collision_report_raw9.txt", 'r')
collisions_raw = [line.split() for line in f.readlines()]
times = []
targ_hashes = []
targ_masses = []
proj_hashes = []
proj_masses = []
mlr_over_mtarg = []
collision_types = [] #gonna need to fix this for head-on smashed collisions - need to decide if they're 2 or 3 based off mlr/mt
frags_raw = []


print(collisions_raw[4])
for lst in collisions_raw:
    if lst[0] == 'TIME':
        times.append(float(lst[3]))
    if lst[0] == 'Target':
        targ_hashes.append(int(lst[4]))
        targ_masses.append(float(lst[5]))
    if lst[0] == 'Projectile':
        proj_hashes.append(int(lst[4]))
        proj_masses.append(float(lst[5]))
    if lst[0] == 'Mlr/Mt:':
        mlr_over_mtarg.append(float(lst[1]))
    if lst[0] == 'COLLISION':
        if lst[2] == 'ELASTIC' and lst[3] == 'BOUNCE':
            collision_types.append(0)
        if lst[2] == 'SIMPLY' and lst[3] == 'MERGED':
            collision_types.append(1)
        if lst[2] == 'EFFECTIVELY' and lst[3] == 'MERGED':
            collision_types.append(1)
        if lst[2] == 'GRAZE' and lst[3] == 'AND' and lst[4] == 'MERGE':
            collision_types.append(1)
        if lst[2] == 'HIT' and lst[3] == 'AND' and lst[4] == 'RUN':
            collision_types.append(2)
        if lst[2] == 'HEAD-ON' and lst[3] == 'SMASHED':
            collision_types.append(3) #Will need to be changed depending
        if lst[2] == 'GRAZING' and lst[3] == 'SMASHED':
            collision_types.append(3) #Will need to be changed depending 
            #might need to add something for super-castasrophic later

bad_chars = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
frag_group = []

for lst in collisions_raw:
    corrected_string = ""
    for i in lst[0]:
        if i not in bad_chars:
            corrected_string+=i
    if corrected_string == 'FRAG':
        frag_group.append(int(lst[3]))
        frag_group.append(float(lst[4]))
    else:
        frags_raw.append(frag_group)
        frag_group=[]

frags = list(filter(None, frags_raw)) #gets rid of empty lists

frag_idx = 0
for i in range(len(collision_types)):
    if collision_types[i] == 2:
        frag_list = frags[frag_idx]
        no_frags = int(len(frag_list)/2)
        frag_masses = [frag_list[(j*2)+1] for j in range(no_frags)]
        tot_frag_mass = sum(frag_masses)
        if proj_masses[i] < tot_frag_mass: #if there's more mass in the fragments than in the projectile, then this is a partial erosion
            collision_types[i] += 1
            print(frag_idx)
            print(times[i])
        frag_idx += 1
    elif collision_types[i] == 3:
        frag_list = frags[frag_idx]
        no_frags = int(len(frag_list)/2)
        frag_masses = [frag_list[(j*2)+1] for j in range(no_frags)]
        tot_frag_mass = sum(frag_masses)
        if proj_masses[i] > tot_frag_mass: #if there's less mass in the fragments than in the projectile, then this is a partial accretion 
            collision_types[i] += -1
        frag_idx += 1
    else:
        continue
        

mlr = []
frag_idx = 0           
for i in range(len(collision_types)):
    if collision_types[i] == 0:
        mlr.append(targ_masses[i])
    elif collision_types[i] == 1:
        mlr.append(targ_masses[i]+proj_masses[i])
    elif collision_types[i] == 2:
        frag_list = frags[frag_idx]
        no_frags = int(len(frag_list)/2)
        frag_masses = [frag_list[(j*2)+1] for j in range(no_frags)]
        tot_frag_mass = sum(frag_masses)
        accreted_mass = abs(proj_masses[i]-tot_frag_mass)
        mlr.append(targ_masses[i]+accreted_mass)
        frag_idx +=1
    elif collision_types[i] == 3 or collision_types[i] == 4:
        frag_list = frags[frag_idx]
        no_frags = int(len(frag_list)/2)
        frag_masses = [frag_list[(j*2)+1] for j in range(no_frags)]
        tot_frag_mass = sum(frag_masses)
        stripped_mass = abs(tot_frag_mass-proj_masses[i])
        mlr.append(targ_masses[i]-stripped_mass)
        frag_idx +=1


frag_idx = 0            
o_f = open("/Users/nferich/Desktop/anna_collision_reports/new_collision_reports/new_collision_report9.txt", 'w')
for i in range(len(times)): #iterates over each index in the individual particle list
    o_f.write(str(times[i]) + ' ')
    o_f.write(str(collision_types[i]) + ' ')
    o_f.write(str(targ_hashes[i]) + ' ')
    o_f.write(str(mlr[i]) + ' ')
    o_f.write(str(proj_hashes[i]) + ' ')
    if collision_types[i] == 2 or collision_types[i] == 3 or collision_types[i] == 4:
        frag_list = frags[frag_idx] #gets list of fragments created in this collision
        for j in range(len(frag_list)):
            o_f.write(str(frag_list[j]) + ' ')
        frag_idx += 1 #moves index up for next collision
    o_f.write('\n')#go to a new line and to a new particle
o_f.close()


            
        
        
