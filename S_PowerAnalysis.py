# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:26:46 2019

@author: Valeria
"""

import iv_save_module as ivs
import iv_utilities_module as ivu
import matplotlib.pyplot as plt

#%% Parameters

this_filename = 'C:\\Users\\Valeria\\OneDrive\\Labo 6 y 7\\Análisis\\Potencia_M_20191018_10\\Resultados.txt'

#%% Load data

this_data, this_header, this_footer = ivs.loadTxt(this_filename)

#%% Plot

# Plot results for the different rods
fig, ax1 = plt.subplots()

# Frequency plot, right axis
ax1.set_xlabel('Repetición')
ax1.set_ylabel('Frecuencia (GHz)', color='tab:red')
ax1.plot(this_data[:,1], 'ro')
ax1.tick_params(axis='y', labelcolor='tab:red')

# Quality factor, left axis
ax2 = ax1.twinx()  # Second axes that shares the same x-axis
ax2.set_ylabel('Tiempo de decaimiento (ps)', color='tab:blue')
ax2.plot(this_data[:,2], 'bx')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
plt.grid(which='both', axis='x')
ax1.tick_params(length=2)
ax1.grid(axis='x', which='both')

# Save plot
ivs.saveFig(this_filename)

#%% Make table

terms_heading = ["Repetición", "F (GHz)", "\u03C4 (ps)", "Q"]
terms_heading = '\t'.join(terms_heading)
terms_table = ['\t'.join([str(element) for element in row]) for row in this_data]
terms_table = '\n'.join(terms_table)
terms_table = '\n'.join([terms_heading, terms_table])
ivu.copy(terms_table)
