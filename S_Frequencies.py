# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 19:51:46 2019

@author: Vall
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import iv_save_module as ivs
import iv_utilities_module as ivu
import iv_analysis_module as iva

#%% PARAMETERS ----------------------------------------------------------------

# Main folder's path
home = r'C:\Users\Luciana\OneDrive\Labo 6 y 7'
# Path to a list of filenames and rods to analize
rods_filename = os.path.join(home, r'Análisis\Rods_LIGO1.txt')
sem_filename = os.path.join(home, r'Muestras\SEM\LIGO1\LIGO1 Geometrías\1\Resultados_LIGO1_1.txt')
desired_frequency = 10 # in GHz
minimum_frequency = 10

#%% LOAD DATA -----------------------------------------------------------------

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
fits_filenames = [ivs.filenameToFitsFilename(file, home) for file in filenames]

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

# Also load data from SEM dimension analysis
other_data, other_header, other_footer = ivs.loadTxt(sem_filename)
other_rods = other_footer['rods']
new_data = []
for r in rods:
    i = other_rods.index(r)
    new_data.append(other_data[i])
sem_data = np.array(new_data)

# Prepare important data for a table 
items = []
for i in range(len(rods)):
    h = '\t'.join(ivu.errorValue(sem_data[i,2], sem_data[i,3]))
    ra = '\t'.join(ivu.errorValue(sem_data[i,4], sem_data[i,5], one_point_scale=True))
    items.append('\t'.join([h, ra, 
                            "{:.2f}".format(fits_data[i,0]), 
                             "{:.1f}".format(fits_data[i,2])]))
del h, ra

# Make OneNote table
heading = '\t'.join(["Longitud (nm)", "Error (nm)", 
                     "Relación de aspecto", "Error",
                     "Frecuencia (GHz)", "Factor de calidad"])
items = ['\t'.join([n, r]) for n, r in zip(rods, items)]
items = '\n'.join(items)
heading = '\t'.join(['Rod', heading])
table = '\n'.join([heading, items])
ivu.copy(table)
del heading, items

# Save all important data to a single file
whole_filename = os.path.join(home, r'Análisis/Resultados_LIGO1_1.txt')
whole_data = np.array([*sem_data[:,:6].T, fits_data[:,0], fits_data[:,2]])
ivs.saveTxt(whole_filename, whole_data, 
            header=["Ancho (nm)", "Error (nm)",
                    "Longitud (nm)", "Error (nm)", 
                    "Relación de aspecto", "Error",
                    "Frecuencia (GHz)", "Factor de calidad"],
            footer=dict(rods=rods),
            overwrite=True)


#%% ANALYSIS ------------------------------------------------------------------

#%% 1) FREQUENCY AND QUALITY FACTOR PER ROD

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
plt.xticks(np.arange(len(rods)), rods, rotation='vertical')
plt.grid(which='both', axis='x')
ax1.tick_params(length=2)
ax1.grid(axis='x', which='both')
#ax1.tick_params(axis='x', labelrotation=90)

#%% 2) FREQUENCY AND LENGTH

#plt.figure()
#plt.plot(1/sem_data[:,2],fits_data[:,0],'o')
#plt.ylabel('Frecuencia (GHz)')
#plt.xlabel('1/Longitud (1/nm)')

length = sem_data[:,2]
young = 45e9 # E = 78 GPa 45
young_2 = 64e9
young_3 = 78e9
density = 19300 # kg/m3
theory = (np.sqrt(young/density) / (2 * length*1e-9)) * 1e-9
theory_2 = (np.sqrt(young_2/density) / (2 * length*1e-9)) * 1e-9
theory_3 = (np.sqrt(young_3/density) / (2 * length*1e-9)) * 1e-9

factor = 0.2 # fraction that is inmersed
omega_0 = np.pi * np.sqrt(young_2/density) / (length*1e-9)
G = 30e9 # Pa
K1 = G * 2.75 # Pa
d = 27e-9 # m           
A = np.pi * (d**2) / 4 # m2
posible_increment = np.sqrt(omega_0**2 + factor*K1/(density*A)) / np.sqrt(omega_0**2)

plt.figure()
plt.plot(length, posible_increment)

factor_4 = 0.2  
theory_4 = np.sqrt(omega_0**2 + factor_4*K1/(density*A)) * 1e-9 / (2 * np.pi)
factor_5 = 0.1
theory_5 = np.sqrt(omega_0**2 + factor_5*K1/(density*A)) * 1e-9 / (2 * np.pi)

plt.figure()
plt.loglog(length, fits_data[:,0],'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.loglog(length, theory, '-', label="{} GPa".format(young/1e9))
plt.loglog(length, theory_2, '-', label="{} GPa".format(young_2/1e9))
plt.loglog(length, theory_3, '-', label="{} GPa".format(young_3/1e9))
plt.loglog(length, theory_4, '-', label="{} en {} GPa".format(factor_4,
                                                              young_3/1e9))
plt.loglog(length, theory_5, '-', label="{} en {} GPa".format(factor_5,
                                                              young_3/1e9))
plt.legend()

rsq, m, b = iva.linearFit(1/length, fits_data[:,0], M = True)
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Inverso de longitud (1/nm)')
young_fit = 4 * density * (m[0]**2)
young_fit_error = 8 * density * m[1] * m[0]
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young_fit, 
                                                        young_fit_error, 
                                                        units="Pa")))

#%% 3) FREQUENCY AND ASPECT RELATION

# Plot results 
fig, ax1 = plt.subplots()

# Frequency vs width plot, lower axis
ax1.set_xlabel('Ancho (nm)', color='tab:red')
ax1.set_ylabel('Frecuencia (GHz)')
ax1.plot(sem_data[:,0], fits_data[:,2], 'ro')
ax1.tick_params(axis='x', labelcolor='tab:red')

# Frequency vs quality factor, upper axis
ax2 = ax1.twiny()  # Second axes that shares the same y-axis
ax2.set_xlabel('Factor de calidad (u.a.)', color='tab:blue')
ax2.plot(sem_data[:,4], fits_data[:,2], 'bx')
ax2.tick_params(axis='x', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
ax1.grid(axis='both')
                     
#%% 4) FREQUENCY AND LENGTH WITHOUT OUTLIERS

# Try it all again, removing outliers
nice_rods = fits_data[:,0] >= minimum_frequency
nice_fits_data = fits_data[nice_rods,:]
nice_sem_data = sem_data[nice_rods, :]

rsq, m, b = iva.linearFit(1/nice_sem_data[:,2], nice_fits_data[:,0])
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Inverso de longitud (1/nm)')
young_fit = 4 * density * (m[0]**2)
young_fit_error = 8 * density * m[1] * m[0]
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young_fit, 
                                                        young_fit_error, 
                                                        units="Pa")))

#%% 5) FREQUENCY AND ASPECT RELATION WITHOUT OUTLIERS

# Plot results 
fig, ax1 = plt.subplots()

# Frequency vs width plot, lower axis
ax1.set_xlabel('Ancho (nm)', color='tab:red')
ax1.set_ylabel('Frecuencia (GHz)')
ax1.plot(nice_sem_data[:,0], nice_fits_data[:,2], 'ro')
ax1.tick_params(axis='x', labelcolor='tab:red')

# Frequency vs quality factor, upper axis
ax2 = ax1.twiny()  # Second axes that shares the same y-axis
ax2.set_xlabel('Factor de calidad (u.a.)', color='tab:blue')
ax2.plot(nice_sem_data[:,4], nice_fits_data[:,2], 'bx')
ax2.tick_params(axis='x', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
ax1.grid(axis='both')

#%% 6) HISTOGRAMS

fig, ax = plt.subplots()

bins_limits = ax.hist(sem_data[:,2])[1]
plt.xlabel("Longitud (nm)")
plt.ylabel("Repeticiones")

#mean_frequencies = []
#for Fi, Ff in zip(bins_limits[:-1], bins_limits[1:]):
#    mean_frequencies = np.mean(fits_data[Fi<=fits_data[:,0]<Ff,0])

fig, ax = plt.subplots()

bins_limits = ax.hist(fits_data[:,0])[1]
plt.xlabel("Frecuencia (GHz)")
plt.ylabel("Repeticiones")

#%% FREQUENCY AND LENGTH FORCED FIT TO SLOPE -1

# Try to fit only the y-interccept, stating that the slope must be -1.
def function(x, b):
    return -x + b
rsq_2, b_2 = iva.nonLinearFit(np.log(sem_data[:,2]), 
                              np.log(fits_data[:,0]), 
                              function)
b_2 = b_2[0]
young_fit_2 = density * (2 * np.exp(b_2[0]))**2
young_fit_error_2 = 2 * density * b_2[1] * (2 * np.exp(b_2[0]))**2 
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young_fit_2, 
                                                        young_fit_error_2, 
                                                        units="Pa")))

#%% 8) BOX PLOTS

plt.figure()
plt.boxplot(fits_data[:,0])

plt.figure()
plt.boxplot(sem_data[:,4])

plt.figure()
plt.boxplot(sem_data[:,2])

#%% SIDE INVESTIGATION :P

Q = np.linspace(12,300,50)
f = 9e9
tau = Q / (np.pi * f)
x = np.linspace(5e9,14e9,500)

curves = np.array([np.imag(f / (f**2 - x**2 - 1j * x / (np.pi * t))) for t in tau])

width = []
for c in curves:
    i = c.argmax()
    x1 = np.argmin(np.abs(c[:i] - max(c) * np.ones(len(c[:i])) / 2))
    x2 = np.argmin(np.abs(c[i:] - max(c) * np.ones(len(c[i:])) / 2)) + i
    width.append(x2-x1)
del c, i, x1, x2
width = np.array(width)

plt.plot(Q, width)
plt.xlabel("Factor de calidad")
plt.ylabel("Ancho de la campana (Hz)")