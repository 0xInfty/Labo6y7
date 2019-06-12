# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 19:51:46 2019

@author: Vall
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import iv_save_module as ivs

# PARAMETERS ------------------------------------------------------------------

# Main folder's path
home = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7'
# Path to a list of filenames and rods to analize
rods_filename = os.path.join(home, 'Análisis\Rods_LIGO1.txt')
desired_frequency = 9.5 # in GHz

# CODE ------------------------------------------------------------------------

# Look for the list of rods and filenames
filenames = [] # Will contain filenames like 'M_20190610_01'
rods = [] # Will contain rods' positions like '1,2'
with open(rods_filename, 'r') as file:
    for line in file:
        if line[0]!='#':
            filenames.append(line.split('\t')[0]) # Save filenames
            rods.append(line.split('\t')[1].split('\n')[0]) # Save rods
    del line

# Now create a list of folders for each filename    
def filenameToFitsFilename(filename):
    """Given a filename 'M_20190610_01', returns path to fits' data"""
    date = filename.split('_')[1] # From 'M_20190610_01' take '20190610'
    date = '-'.join([date[:4], date[4:6], date[6:]]) # Transfrom to '2019-06-10'
    fits_filename = os.path.join(home, 'Mediciones', date, 
                                 'Ajustes', filename+'.txt')
    return fits_filename
fits_filenames = [filenameToFitsFilename(file) for file in filenames]

# Load data from each fit
fits_data = []
fits_footer = []
for file in fits_filenames:
    data, fits_header, footer = ivs.loadTxt(file)
    fits_data.append(data)
    fits_footer.append(footer)
del file, data, footer

# Keep only the fit term that has the closest frequency to the desired one
index = []
for rod, fit in zip(rods, fits_data):
    i = np.argmin(abs(fit[:,0] - desired_frequency*np.ones(fit.shape[0])))
    index.append(i)
del rod, fit, i
fits_data = np.array([fit[i,:] for fit, i in zip(fits_data,index)])
del index

# Plot frequency, characteristic time and quality factor
#plt.figure()
#plt.subplot(3,1,1)
#plt.plot(fits_data[:,0], 'o')
#plt.subplot(3,1,2)
#plt.plot(fits_data[:,2], 'o')
#plt.subplot(3,1,3)
#plt.plot(fits_data[:,1], 'o')

# Plot results for the different rods
fig, ax1 = plt.subplots()

# Frequency plot, right axis
ax1.set_xlabel('Antena')
ax1.set_ylabel('Frecuencia (GHz)', color='tab:red')
ax1.plot(fits_data[:,0], 'ro')
ax1.tick_params(axis='y', labelcolor='tab:red')

# Quality factor, left axis
ax2 = ax1.twinx()  # Second axes that shares the same x-axis
ax2.set_ylabel('Factor de calidad (u.a.)', color='tab:blue')
ax2.plot(fits_data[:,2], 'bx')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
plt.xticks(np.arange(len(rods)), rods)
plt.grid(which='both', axis='x')
ax1.tick_params(length=2)
#for tick in ax1.axes.get_xticklabels():
#tick.set_visible(False)
ax1.grid(axis='x', which='both')
ax1.tick_params(axis='x', labelrotation=-90)

other_data, other_header, other_footer = ivs.loadTxt(r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7\Muestras\SEM\LIGO1\LIGO1 Geometrías\1\Resultados_LIGO1_1.txt')
other_rods = other_footer['rods']
new_data = []
for r in rods:
    i = other_rods.index(r)
    new_data.append(other_data[i])
    print(i)
sem_data = np.array(new_data)

#%% ANALYSIS WITH ANDREA ------------------------------------------------------

#plt.figure()
#plt.plot(1/sem_data[:,2],fits_data[:,0],'o')
#plt.ylabel('Frecuencia (GHz)')
#plt.xlabel('1/Longitud (1/nm)')

length = sem_data[:,2]
young = 42e9 # E = 78 GPa
young_2 = 79e9
density = 19300 # kg/m3
theory = (np.sqrt(young/density) / (2 * length*1e-9)) * 1e-9
theory_2 = (np.sqrt(young_2/density) / (2 * length*1e-9)) * 1e-9

plt.figure()
plt.loglog(sem_data[:,2],fits_data[:,0],'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.loglog(length, theory, 'x')
plt.loglog(length, theory_2, 'x')