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
home = r'C:\Users\Vall\OneDrive\Labo 6 y 7'
desired_frequency = 10 # in GHz
minimum_frequency = 10

# Path to a list of filenames and rods to analize
rods_filename = os.path.join(home, r'Análisis\Rods_LIGO1.txt')
sem_series = ['LIGO1_1']
sem_short_series = lambda series : '{}'#series.split('_')[1]+' {}'
name = 'LIGO1'

"""
rods_filename = os.path.join(home, r'Análisis\Rods_M135.txt')
sem_series = ['M135_5_1D', 'M135_7B_1D']
sem_short_series = lambda series : series.split('_')[1]+' {}'
name = 'M135'
"""

def filenameToSEMFilename(series, home=home):
    
    """Given a series 'M135_7B_1D', returns path to SEM data"""
    
    filename = 'Resultados_SEM_{}.txt'.format(series)
    series = series.split('_') # From 'M_20190610_01' take '20190610'
    sem_filename = os.path.join(home, 'Muestras\SEM', *series, filename)
    
    return sem_filename

sem_filename = [filenameToSEMFilename(s) for s in sem_series]

#%% LOAD DATA -----------------------------------------------------------------

# Look for the list of rods and filenames
filenames = [] # Will contain filenames like 'M_20190610_01'
rods = [] # Will contain rods' positions like '1,2'
with open(rods_filename, 'r') as file:
    for line in file:
        if line[0]!='#':
            filenames.append(line.split('\t')[0]) # Save filenames
            aux = line.split('\t')[1:]
            aux = r' '.join(aux)
            rods.append(aux.split('\n')[0]) # Save rods
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
fits_new_data = []
for rod, fit in zip(rods, fits_data):
    try:
        i = np.argmin(abs(fit[:,0] - desired_frequency*np.ones(fit.shape[0])))
        fits_new_data.append([*fit[i,:]])
    except IndexError:
        fits_new_data.append([*fit])
fits_data = np.array(fits_new_data)
del rod, fit, i, fits_new_data

# Also load data from SEM dimension analysis
sem_data = []
other_rods = []
for sf, s in zip(sem_filename, sem_series):
    d, other_header, f = ivs.loadTxt(sf)
    r = [sem_short_series(s).format(rs) for rs in f['rods']]
    sem_data.append(d)
    other_rods = other_rods + r
del d, f, r

other_data = []
for s in sem_data:
    for si in s:
        other_data.append([*si])
other_data =  np.array(other_data)
del sem_data

sem_data = []
for r in rods:
    i = other_rods.index(r)
    sem_data.append(other_data[i])
sem_data = np.array(sem_data)

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
whole_filename = os.path.join(home, r'Análisis/Resultados_Totales_{}.txt'.format(name))
whole_data = np.array([*sem_data[:,:6].T, fits_data[:,0], fits_data[:,2]])
ivs.saveTxt(whole_filename, whole_data.T, 
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
ax1.tick_params(axis='x', labelrotation=90)

#%% 2A) FREQUENCY AND LENGTH
# --> Try out some known values

#plt.figure()
#plt.plot(1/sem_data[:,2],fits_data[:,0],'o')
#plt.ylabel('Frecuencia (GHz)')
#plt.xlabel('1/Longitud (1/nm)')

# Parameters
density = 19300 # kg/m3
G = 30e9 # Pa
K1 = G * 2.75 # Pa
diameter = 27e-9 # m
# Why is this only lateral surface?

# Data
length = sem_data[:,2]
young = [45e9, 64e9, 78e9]
factor = [0, .1, .2, .3, .4] # surface fraction that is inmerse

# Theory models
def area(length, diameter=diameter):
    A = np.pi * (diameter**2) / 4#  + np.pi * diameter * length # m2
    # Why is this only lateral surface?
    return A

def f_simple(length, young):
    f_0 = (np.sqrt(young/density) / (2 * length*1e-9)) * 1e-9
    return f_0

def f_complex(length, young, factor, density=density, K1=K1):
    f_0 = f_simple(length, young)
    f = np.sqrt( (2*np.pi*f_0*1e9)**2 + factor*K1/(density*area(length)) ) 
    f = f * 1e-9/(2*np.pi)
    return f

# Theory predictions
f_0 = np.array([f_simple(l, y) for l in length for y in young])
f_0 = f_0.reshape([len(length), len(young)])
f = np.array([f_complex(l, young[1], c) for l in length for c in factor])
f = f.reshape([len(length), len(factor)])

#posible_increment = np.sqrt(omega_0**2 + factor*K1/(density*A)) / np.sqrt(omega_0**2)
#plt.figure()
#plt.plot(length, posible_increment)

plt.figure()
plt.title('Modelo simple')
plt.loglog(length, fits_data[:,0],'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in f_0.T: plt.loglog(length, freq, '-')
plt.legend(["Datos"] + ["{} GPa".format(y/1e9) for y in young])

plt.figure()
plt.title('Modelo complejo con Young {} GPa'.format(young[1]/1e9))
plt.loglog(length, fits_data[:,0],'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in f.T: plt.loglog(length, freq, '-')
plt.legend(["Datos"] + ["Sumergido {:.0f}%".format(c*100) for c in factor])

#%% 2B) FREQUENCY AND LENGTH
# --> Try a linear fit

rsq, m, b = iva.linearFit(1/length, fits_data[:,0], M = True)
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Inverso de longitud (1/nm)')
young_fit = 4 * density * (m[0]**2)
young_fit_error = np.abs(8 * density * m[1] * m[0])
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young_fit, 
                                                        young_fit_error, 
                                                        units="Pa")))

#%% 2C) FREQUENCY AND LENGTH
# --> Try a linear fit with a forced slope -1

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

#%% 2D) FREQUENCY AND LENGTH
# --> Try a nonlinear fit with complex model

# Parameters
density = 19300 # kg/m3
G = 30e9 # Pa
K1 = G * 2.75 # Pa
diameter = 27e-9 # m
# Why is this only lateral surface?

# Data
young = [45e9, 64e9, 78e9]
factor = [0, .1, .2, .3, .4] # surface fraction that is inmerse

# Theory models
def simple_function(length, young):
    f_0 = (np.sqrt(young/density) / (2 * length*1e-9)) * 1e-9
    return f_0    

def complex_function(length, young, factor):
    A = np.pi * (diameter**2) / 4#  + np.pi * diameter * length # m2
    f_0 = (np.sqrt(young/density) / (2 * length*1e-9)) * 1e-9
    f = np.sqrt( (2*np.pi*f_0*1e9)**2 + factor*K1/(density * A) ) 
    f = f * 1e-9/(2*np.pi)
    return f

# Mohsen initial values (64e9, .3)
rsq_nl_c, parameters_nl_c = iva.nonLinearFit(
        sem_data[:,2], 
        fits_data[:,0],
        complex_function)

#%% 3) FREQUENCY VS Q AND WIDTH

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

#%% 4) HISTOGRAMS

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


#%% 5) BOX PLOTS

plt.figure()
plt.boxplot(fits_data[:,0])

plt.figure()
plt.boxplot(sem_data[:,4])

plt.figure()
plt.boxplot(sem_data[:,2])

#%% *) RAMAN-LIKE SPECTRUM SIDE INVESTIGATION

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