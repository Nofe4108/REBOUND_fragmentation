"""
Produces Figures 11 and 12 from Ferich et al. (IN PREP).
The plots focus on planets produced by each REBOUND simulation,
which are objects that have a mass greater than 0.093 Mearth
"""

import numpy as np
import time
import matplotlib.pyplot as plt
from astropy import units as u
from Differentiated_Body_Composition_Tracker_for_Paper import track_composition

start_time = time.time() #Timer to see how long running this code takes

# Constants
mass_conversion = ((1.0*u.Msun).to(u.Mearth)).value # Number that converts Msun to Mearth

# File Pathways
collision_file_pw = "new_collision_reports/new_collision_report"
composition_input_file_pw = "DBCT_input/uni_DBCT_input"
impact_parameter_file_pw = "impact_parameters/impact_parameters"
ejection_file_pw = "ejections/ejections"
output_file_pw = "DBCT_output/uni_DBCT_output"

file_range = np.arange(1,51,1) # Number of simulations the code produces data for

########### DATA COLLECTION AND PLOTTING FOR FIGURE 1 ##################

all_final_cmfs = []
final_average_cmfs = []
min_final_cmfs = []
max_final_cmfs = []

min_ejecta_cmfs = np.arange(0, 1.1, 0.1) #mimimum fraction of core material in ejecta
max_ejecta_cmfs = np.arange(1, 2, 1) #maximum fraction of core material in ejecta

# This loop runs the DBCT with different combinations of fmin and fmax and extracts CMF data for planets
for max_ejecta_cmf in max_ejecta_cmfs:
    for min_ejecta_cmf in min_ejecta_cmfs:
        final_cmfs = []
        for i in file_range:
            collision_file = collision_file_pw + str(i) + ".txt"
            composition_input_file = composition_input_file_pw + ".txt"
            impact_parameter_file = impact_parameter_file_pw + str(i) + ".txt"
            ejection_file = ejection_file_pw + str(i) + ".txt"
            final_composition = track_composition(collision_file, composition_input_file, ejection_file, impact_parameter_file, min_ejecta_cmf, max_ejecta_cmf)
            for obj in final_composition:
                if obj[1]*mass_conversion > 0.093: #if the object mass is bigger than original embryo mass
                    final_cmfs.append(obj[2])
        all_final_cmfs.append(final_cmfs)
        final_average_cmfs.append(sum(final_cmfs)/len(final_cmfs))
        min_final_cmfs.append(min(final_cmfs))
        max_final_cmfs.append(max(final_cmfs))

# Removes the CMF of planet with the most extreme CMF to prevent skewing
for i in range(len(all_final_cmfs)):
    all_final_cmfs[i].pop(26)

# Figure 11      
fig1, ax1 = plt.subplots(figsize=(6,5))
labels = [str(round(num, 1)) for num in np.arange(0.0, 1.1, 0.1)]
ax1.axhline(0.3, label = 'Initial CMF', color = 'tab:orange', linestyle = '-', alpha=0.6)
fig1 = plt.violinplot(all_final_cmfs, showmedians=True, widths=0.7)
ax1.set_xlabel('$f_{min}$', fontsize='x-large')
ax1.set_ylabel('Final Planet CMF', fontsize='large')
ax1.set_title('$f_{max}$ = 1.0', fontsize='x-large')
ax1.grid(color='grey', linestyle='-', axis='y', alpha=0.7)
ax1.minorticks_on()
ax1.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
ax1.set_yticks(np.arange(0.00, 0.85, 0.0125)) #.4
ax1.set_ylim([0.225, 0.375])
ax1.legend(loc='lower right')
plt.savefig("graphs/varying_fmin_violin_plot.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/varying_fmin_violin_plot.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/varying_fmin_violin_plot.png', bbox_inches='tight', dpi=300)

########### DATA COLLECTION AND PLOTTING FOR FIGURE 2 ##################

all_final_cmfs = []
final_average_cmfs = []
min_final_cmfs = []
max_final_cmfs = []

min_ejecta_cmfs = np.arange(0, 1, 1) #mimimum fraction of core material in ejecta
max_ejecta_cmfs = np.arange(0, 1.1, 0.1) #maximum fraction of core material in ejecta

# This loop runs the DBCT with different combinations of fmin and fmax and extracts CMF data for planets
for min_ejecta_cmf in min_ejecta_cmfs:
    for max_ejecta_cmf in max_ejecta_cmfs:
        final_cmfs = []
        for i in file_range:
            collision_file = collision_file_pw + str(i) + ".txt"
            composition_input_file = composition_input_file_pw + ".txt"
            impact_parameter_file = impact_parameter_file_pw + str(i) + ".txt"
            ejection_file = ejection_file_pw + str(i) + ".txt"
            final_composition = track_composition(collision_file, composition_input_file, ejection_file, impact_parameter_file, min_ejecta_cmf, max_ejecta_cmf)
            for obj in final_composition:
                if obj[1]*mass_conversion > 0.093: #if the object mass is bigger than original embryo mass
                    final_cmfs.append(obj[2])
        all_final_cmfs.append(final_cmfs)
        final_average_cmfs.append(sum(final_cmfs)/len(final_cmfs))
        min_final_cmfs.append(min(final_cmfs))
        max_final_cmfs.append(max(final_cmfs))

# Removes the CMF of planet with the most extreme CMF to prevent skewing
for i in range(len(all_final_cmfs)):
    all_final_cmfs[i].pop(26)

# Figure 11      
fig1, ax1 = plt.subplots(figsize=(6,5))
labels = [str(round(num, 1)) for num in np.arange(0.0, 1.1, 0.1)]
ax1.axhline(0.3, label = 'Initial CMF', color = 'tab:orange', linestyle = '-', alpha=0.6)
fig1 = plt.violinplot(all_final_cmfs, showmedians=True, widths=0.7)
ax1.set_xlabel('$f_{max}$', fontsize='x-large')
ax1.set_ylabel('Final Planet CMF', fontsize='large')
ax1.set_title('$f_{min}$ = 0.0', fontsize='x-large')
ax1.grid(color='grey', linestyle='-', axis='y', alpha=0.7)
ax1.minorticks_on()
ax1.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
ax1.set_yticks(np.arange(0.00, 0.85, 0.0125)) #.4
ax1.set_ylim([0.225, 0.375])
ax1.legend(loc='lower right')
plt.savefig("graphs/varying_fmax_violin_plot.pdf", bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/varying_fmax_violin_plot.eps', bbox_inches='tight', pad_inches=0.01)
plt.savefig('graphs/varying_fmax_violin_plot.png', bbox_inches='tight', dpi=300)

print("--- %s seconds ---" % (time.time() - start_time))