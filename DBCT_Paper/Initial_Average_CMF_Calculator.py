# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 18:23:10 2024

@author: nferi

Calculates the initial average CMF for a particular
distribution
"""
import numpy as np
import matplotlib.pyplot as plt
from Differentiated_Body_Composition_Tracker_for_Paper import organize_compositions

final_cmfs = []

composition_input_file = "DBCT_input/reversed_3step_DBCT_input.txt"

f = open(composition_input_file, 'r')
init_compositions = [line.split() for line in f.readlines()] # Reads each line in file and splits its value into an array (hash - mass - composition fractions) (all values are strings)
compositions = organize_compositions(init_compositions) # Organzies the data from the file
init_cmfs = [obj[2] for obj in compositions]
init_avg_cmf = np.average(init_cmfs)

print("Initial Average CMF:", str(init_avg_cmf))
