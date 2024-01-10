# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:51:08 2023

@author: nferi
"""
"""This code will produce the final hashes, masses, and CMFs for collisional data
from Childs et al. (2022). Different input files can be entered for the initial
objects in the disc"""

import numpy as np
import time
from Differentiated_Body_Composition_Tracker_for_Paper import track_composition
from Differentiated_Body_Composition_Tracker_for_Paper import write_output

start_time = time.time() #Timer to see how long running this code takes

min_ejecta_cmf = 0.0 #mimimum fraction of core material in ejecta
max_ejecta_cmf = 1.0 #maximum fraction of core material in ejecta

collision_file_pw = "new_collision_reports/new_collision_report"
composition_input_file_pw = "DBCT_input/uni_DBCT_input"
impact_parameter_file_pw = "impact_parameters/impact_parameters"
ejection_file_pw = "ejections/ejections"
output_file_pw = "DBCT_output/uni_DBCT_output"

file_range = np.arange(1,51,1) # Number of simulations the code produces data for

for i in file_range:
    print(i)
    collision_file = collision_file_pw + str(i) + ".txt"
    composition_input_file = composition_input_file_pw + ".txt"
    impact_parameter_file = impact_parameter_file_pw + str(i) + ".txt"
    ejection_file = ejection_file_pw + str(i) + ".txt"
    output_file = output_file_pw + str(i) + ".txt"

    write_output(track_composition(collision_file, composition_input_file, ejection_file, impact_parameter_file, min_ejecta_cmf, max_ejecta_cmf), output_file)
    
    
    
    