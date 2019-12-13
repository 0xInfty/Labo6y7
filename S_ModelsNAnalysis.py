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
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
desired_frequency = 8 # in GHz

# Path to a list of filenames and rods to analize

rods_filename = os.path.join(home, r'Análisis\Rods_LIGO5bis.txt')
sem_series = ['LIGO5bis_1']
sem_short_series = lambda series : '{}'#series.split('_')[1]+' {}'
name = 'LIGO5bis'
"""
rods_filename = os.path.join(home, r'Análisis\Rods_M135.txt')
sem_series = ['M135_5_1D', 'M135_7B_1D']
sem_short_series = lambda series : series.split('_')[1]+' {}'
name = 'M135'
"""
# Some function to manege filenames
def filenameToSEMFilename(series, home=home):
    
    """Given a series 'M135_7B_1D', returns path to SEM data"""
    
    filename = 'Resultados_SEM_{}.txt'.format(series)
    series = series.split('_') # From 'M_20190610_01' take '20190610'
    sem_filename = os.path.join(home, 'Muestras\SEM', *series, filename)
    
    return sem_filename

def paramsFilename(series, home=home):
    
    """Given a series name like 'LIGO1', returns path to parameters"""
    
    filename = os.path.join(home, r'Análisis/Params_{}.txt'.format(series))
    
    return filename

def figsFilename(fig_name, series='', home=home):
    
    """Given a fig_name 'DifCuadráticaMedia', returns path to fig"""
    
    if series!='':
        series = '_{}'.format(series)
    base = os.path.join(home, r'Análisis/ModelsNAnalysis'+series)
    if not os.path.isdir(base):
        os.makedirs(base)
    
    filename = os.path.join(base, fig_name+'.png')
    
    return filename

#%%

# Physics' Parameters
density = 19.3e3 # kg/m3 for gold
Shear = np.mean([30.8e9, 32.3e9]) # Pa for fused silica
diameter = 27.7e-9 # m for rods
midlength = 85e-9 # m for rods
Viscosity = 2e-3 # Pa/s for gold
Young = np.mean([71.2e9, 74.8e9])  # Pa for fused silica
density_s = np.mean([2.17e3, 2.22e3]) # kg/m3 for fused silica
area = np.pi * diameter**2 / 4
K1 = Shear * 2.75 # Pa
K2 = np.pi * diameter * np.sqrt(density_s * Shear) # Pa.s (viscosity's units)

# Space to save results
young = {}
#k1 = {}
#k2 = {}
factor = {}
chi_squared = {}

# Theory models
def f_simple(length, young):
    f_0 = (np.sqrt(young/density) / (2 * length))
    return f_0

def f_mid(length, young):
    f_0 = f_simple(length, young)
    beta = ( Viscosity / (length * density) )**2 / 2
    f = np.sqrt(f_0**2 - beta**2 / 4)
    return f

def f_full(length, young):
    f_0 = f_simple(length, young)
    beta = ( Viscosity / (length * density) )**2 / 2
    K1_term = K1 / ( np.pi**2 * density * area )
    K2_subterm = K2 / ( 2 * np.pi * density * area )
    f = np.sqrt(f_0**2 + K1_term/4 - (K2_subterm + beta/np.pi)**2/4 )
    return f

def f_andrea(length, young, factor):
    f_0 = f_simple(length, young)
    f = np.sqrt( (2*np.pi*f_0)**2 + factor*K1/(density*area) ) 
    f = f /(2*np.pi)
    return f

def f_iv(length, young, factor):
    f_0 = f_simple(length, young)
    K1_term = factor * K1 / ( np.pi**2 * density * area )
    K2_subterm = factor * K2 / ( 4 * np.pi * density * area )
    f = np.sqrt(f_0**2 + K1_term/4 - (K2_subterm)**2 )
    return f

def f_test(length, young, shear):
    K1 = 2.75 * shear
    K2 = np.pi * diameter * np.sqrt(density_s * shear)
    f_0 = f_simple(length, young)
    beta = ( Viscosity / (length * density) )**2 / 2
    K1_term = K1 / ( np.pi**2 * density * area )
    K2_subterm = K2 / ( 2 * np.pi * density * area )
    f = np.sqrt(f_0**2 + K1_term/4 - (K2_subterm + beta/np.pi)**2/4 )
    return f

#def f_free(length, young, K1, K2):
#    f_0 = f_simple(length, young)
#    beta = ( Viscosity / (length * density) )**2 / 2
#    K1_term = K1 / ( np.pi**2 * density * area )
#    K2_subterm = K2 / ( 2 * np.pi * density * area )
#    f = np.sqrt(f_0**2 + K1_term/4 - (K2_subterm + beta/np.pi)**2/4 )
#    return f

def tau_simple(length, viscosity):
    tau = 2 * (length * density / (np.pi * viscosity))**2
    return tau

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
            del aux
    del line
del rods_filename

# Then load parameters
params_filenames = [] # Will contain filenames like 'M_20190610_01'
params = {} # Will contain parameters
amplitude = []
power = []
wavelength = []
spectral_width = []
with open(paramsFilename(name), 'r') as file:
    for line in file:
        if line[0]!='#':
            params_filenames.append(line.split('\t')[0])
            amplitude.append(float(line.split('\t')[1]))
            power.append(float(line.split('\t')[2]))
            wavelength.append(float(line.split('\t')[3]))
            spectral_width.append(float(line.split('\t')[-1]))
    del line
params = np.array([amplitude, power, wavelength, spectral_width]).T
index = [params_filenames.index(f) for f in filenames]
params = params[index,:]
params_header = ['Amplitud (mVpp)', 'Potencia Pump post-MOA (muW)', 
                 'Longitud de onda (nm)', 'Ancho medio de la campana (nm)']
del params_filenames, index, amplitude, power, wavelength, spectral_width

# Now create a list of folders for each filename    
fits_filenames = [ivs.filenameToFitsFilename(file, home) for file in filenames]

# Load data from each fit
fits_data = []
fits_footer = []
for file in fits_filenames:
    data, fits_header, footer = ivs.loadTxt(file)
    fits_data.append(data)
    fits_footer.append(footer)
del file, data, footer, fits_filenames

# Keep only the fit term that has the closest frequency to the desired one
fits_new_data = []
for rod, fit in zip(rods, fits_data):
    try:
        i = np.argmin(abs(fit[:,0] - desired_frequency*np.ones(fit.shape[0])))
        fits_new_data.append([*fit[i,:]])
    except IndexError:
        fits_new_data.append([*fit])
fits_data = np.array(fits_new_data)
frequency = fits_data[:,0]*1e9 # Hz
damping_time = fits_data[:,1]*1e-12 # s
quality_factor = fits_data[:,2]
del rod, fit, i, fits_new_data

# Also create a list of folders for SEM filenames
sem_filenames = [filenameToSEMFilename(s) for s in sem_series]

# Then load data from SEM dimension analysis
sem_data = []
sem_footers = []
sem_rods = []
for sf, s in zip(sem_filenames, sem_series):
    d, sem_header, f = ivs.loadTxt(sf)
    r = [sem_short_series(s).format(rs) for rs in f['rods']]
    sem_data.append(d)
    sem_footers = sem_footers + [f]
    sem_rods = sem_rods + r
del sf, s, d, f, r, sem_filenames

other_data = []
for s in sem_data:
    for si in s:
        other_data.append([*si])
other_data =  np.array(other_data)
del s, si

index = [sem_rods.index(r) for r in rods]
sem_data = [other_data[i] for i in index]
sem_data = np.array(sem_data)
length = sem_data[:,2] * 1e-9 # m
width = sem_data[:,0] * 1e-9 # m
del other_data, index, sem_rods

# Now we can filter the results
index = np.argsort(frequency) # Remove the two lowest frequencies
length = length[index[2:]]
width = width[index[2:]]
frequency = frequency[index[2:]]
damping_time = damping_time[index[2:]]
quality_factor = quality_factor[index[2:]]
del index

# Since I'll be analysing frequency vs length mostly...
index = np.argsort(length)
length = length[index]
width = width[index]
frequency = frequency[index]
damping_time = damping_time[index]
quality_factor = quality_factor[index]
del index

#%%

# Prepare important data for a table 
items = []
for i in range(len(rods)):
    h = '\t'.join(ivu.errorValue(sem_data[i,2], sem_data[i,3]))
    ra = '\t'.join(ivu.errorValue(sem_data[i,4], sem_data[i,5], one_point_scale=True))
    items.append('\t'.join([h, ra, 
                            "{:.2f}".format(fits_data[i,0]), 
                             "{:.1f}".format(fits_data[i,2])]))
del i, h, ra

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
whole_data = np.array([*sem_data[:,:6].T, fits_data[:,0], 
                       fits_data[:,1], fits_data[:,2]])
ivs.saveTxt(whole_filename, whole_data.T, 
            header=["Ancho (nm)", "Error (nm)",
                    "Longitud (nm)", "Error (nm)", 
                    "Relación de aspecto", "Error",
                    "Frecuencia (GHz)", "Tiempo de decaimiento (ps)",
                    "Factor de calidad"],
            footer=dict(rods=rods, filenames=filenames),
            overwrite=True)
del whole_data, whole_filename

#%% ANALYSIS ------------------------------------------------------------------

#%% 1A) FREQUENCY AND QUALITY FACTOR PER ROD

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

# Save plot
plt.savefig(figsFilename('FyQvsRod', name), bbox_inches='tight')

#%% 1B) FREQUENCY AND LENGTH PER ROD

# Plot results for the different rods
fig, ax1 = plt.subplots()

# Frequency plot, right axis
ax1.set_xlabel('Antena')
ax1.set_ylabel('Frecuencia (GHz)', color='tab:red')
ax1.plot(fits_data[:,0], 'ro')
ax1.tick_params(axis='y', labelcolor='tab:red')

# Quality factor, left axis
ax2 = ax1.twinx()  # Second axes that shares the same x-axis
ax2.set_ylabel('Longitud (nm)', color='tab:blue')
ax2.plot(sem_data[:,2], 'bx')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
plt.xticks(np.arange(len(rods)), rods, rotation='vertical')
plt.grid(which='both', axis='x')
ax1.tick_params(length=2)
ax1.grid(axis='x', which='both')
ax1.tick_params(axis='x', labelrotation=90)

# Save plot
plt.savefig(figsFilename('FyLvsRod', name), bbox_inches='tight')

#%% 1C) LENGTH AND WIDTH PER ROD

# Plot results for the different rods
fig, ax1 = plt.subplots()

# Frequency plot, right axis
ax1.set_xlabel('Antena')
ax1.set_ylabel('Longitud (nm)', color='tab:red')
ax1.plot(sem_data[:,2], 'ro')
ax1.tick_params(axis='y', labelcolor='tab:red')

# Quality factor, left axis
ax2 = ax1.twinx()  # Second axes that shares the same x-axis
ax2.set_ylabel('Ancho (nm)', color='tab:blue')
ax2.plot(sem_data[:,0], 'bx')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
plt.xticks(np.arange(len(rods)), rods, rotation='vertical')
plt.grid(which='both', axis='x')
ax1.tick_params(length=2)
ax1.grid(axis='x', which='both')
ax1.tick_params(axis='x', labelrotation=90)

# Save plot
plt.savefig(figsFilename('LyAvsRod', name), bbox_inches='tight')

#%% 2A) FREQUENCY AND LENGTH
# --> Try out some known values

# Data
young_predict = [0, 64e9, 78e9, 42e9]
young_predict_select = 64e9
factor_predict = [0, .1, .2, 1] # fraction that represents bound
shear_predict = [10e9, Shear, 40e9, Young]

# Theory predictions
freq_simple = np.array([f_simple(length, y) for y in young_predict]).T
freq_full = np.array([f_full(length, y) for y in young_predict]).T
freq_andrea = np.array([f_andrea(length, young_predict_select, c) 
                for c in factor_predict]).T
freq_iv = np.array([f_iv(length, young_predict_select, c) 
                 for c in factor_predict]).T
freq_test = np.array([f_test(length, 0, s) for s in shear_predict]).T

# Make a plot for the simpler model
plt.figure()
ax = plt.subplot()
plt.title('Modelo simple')
plt.loglog(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in freq_simple.T: plt.loglog(length*1e9, freq*1e-9, '-')
del freq
plt.legend(["Datos"] + ["{} GPa".format(y/1e9) for y in young_predict])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('Simple_Predict', name), bbox_inches='tight')

# Make a plot for the complex model
plt.figure()
ax = plt.subplot()
plt.title('Modelo completo')
plt.loglog(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in freq_full.T: plt.loglog(length*1e9, freq*1e-9, '-')
del freq
plt.legend(["Datos"] + ["{} GPa".format(y/1e9) for y in young_predict])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('Full_Predict', name), bbox_inches='tight')

# Make one too for the Andrea model
plt.figure()
ax = plt.subplot()
plt.title('Modelo completo aproximado con factor con Young {} GPa'.format(young_predict_select/1e9))
plt.loglog(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in freq_andrea.T: plt.loglog(length*1e9, freq*1e-9, '-')
del freq
plt.legend(["Datos"] + ["Factor {:.0f}%".format(c*100) for c in factor_predict])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('Andrea_{}_GPa'.format(young_predict_select/1e9), 
                         name), 
            bbox_inches='tight')

# Make one too for the IV model
plt.figure()
ax = plt.subplot()
plt.title('Modelo completo con factor con Young {} GPa'.format(young_predict_select/1e9))
plt.loglog(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in freq_iv.T: plt.loglog(length*1e9, freq*1e-9, '-')
del freq
plt.legend(["Datos"] + ["Factor {:.0f}%".format(c*100) for c in factor_predict])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('IV_{}_GPa'.format(young_predict_select/1e9), 
                         name), 
            bbox_inches='tight')

# Make one too for the test model
plt.figure()
ax = plt.subplot()
plt.title('Modelo test Shear con Young 0 GPa')
plt.loglog(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in freq_test.T: plt.loglog(length*1e9, freq*1e-9, '-')
del freq
plt.legend(["Datos"] + ["Shear {:.0f} GPa".format(s/1e9) for s in shear_predict])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('Full_0_GPa', name), 
            bbox_inches='tight')

"""LIGO1: Decidimos que esto no es necesario si hacemos ajustes"""

del young_predict, factor_predict, young_predict_select, shear_predict
del freq_simple, freq_full, freq_andrea, freq_iv, freq_test

#%% 2B) FREQUENCY AND LENGTH
# --> Try a linear fit on f vs 1/L, which corresponds to the simple model

#rsq, m, b = iva.linearFit(np.log(length), np.log(frequency), showplot=False)
#plt.figure()
#plt.plot(np.log(length), np.log(fits_data[:,0]), '.')
#plt.plot(np.log(length), m[0]*np.log(length)+b[0], 'r-')
#plt.ylabel('Frecuencia (GHz)')
#plt.xlabel('Inverso de longitud (1/nm)')
#young_fit = 4 * density * (m[0]**2)
#young_fit_error = np.abs(8 * density * m[1] * m[0])
#print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young_fit, 
#                                                        young_fit_error, 
#                                                        units="Pa")))

"""LIGO 1: Da pendiente -1.5 aproximadamente"""

#%% 2C) FREQUENCY AND LENGTH
# --> Try a linear fit with a forced slope -1, even closer to simple model

#def f_simple_loglog(loglength, young):
#    return -loglength + np.log( np.sqrt(young/density) / 2 )
#
#rsq, young = iva.nonLinearFit(np.log(length), 
#                              np.log(frequency), 
#                              f_simple_loglog, showplot=False)
#young = young[0]
#print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young[0], 
#                                                        young[1], 
#                                                        units="Pa")))

"""LIGO1: Esto funciona, pero no hace falta hacerlo tan rebuscado porque puedo 
hacer un ajuste no lineal por cuadrados mínimos directamente de la función."""

#%% 2D) FREQUENCY AND LENGTH
# --> Try a nonlinear fit directly using the simple model

young['simple'] = iva.nonLinearFit(length, frequency, 
                                   f_simple, showplot=False)[-1][0]
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(
            young['simple'][0], 
            young['simple'][1], 
            units="Pa")))
chi_squared['simple'] = sum( (f_mid(length, young['simple'][0]) - frequency)**2 ) 
chi_squared['simple'] = chi_squared['simple'] / len(length)

plt.figure()
ax = plt.subplot()
plt.title('Ajuste modelo simple')
plt.plot(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.plot(length*1e9, 1e-9*f_simple(length, young['simple'][0]), '-r')
plt.legend(["Datos","Ajuste"])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('Simple_Fit', name), bbox_inches='tight')

#%% 2E) FREQUENCY AND LENGTH
# --> Try a nonlinear fit directly using the mid model

#young['mid'] = iva.nonLinearFit(length, frequency, 
#                                f_mid, showplot=False)[-1][0]
#print(r"Módulo de Young: {}".format(ivu.errorValueLatex(
#            young['mid'][0], 
#            young['mid'][1], 
#            units="Pa")))
#chi_squared['mid'] = sum( (f_mid(length, young['mid'][0]) - frequency)**2 ) 
#chi_squared['mid'] = chi_squared['mid'] / len(length)

"""LIGO1: Es muy similar al modelo simple, por lo que el término beta no 
influye demasiado"""

#%% 2F) FREQUENCY AND LENGTH
# --> Try a nonlinear fit directly using the full model

young['full'] = {}
chi_squared['full'] = {}
 
rsq, y = iva.nonLinearFit(length, frequency, f_full, 
                          bounds=([0], [np.infty]), 
                          showplot=False)
young['full'] = y[0]
del y
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young['full'][0], 
                                                        young['full'][1], 
                                                        units="Pa")))

chi_squared['full'] = sum( (f_full(length, young['full'][0]) - frequency)**2 ) 
chi_squared['full'] = chi_squared['full'] / len(length)

plt.figure()
ax = plt.subplot()
plt.title('Ajuste modelo completo')
plt.plot(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.plot(length*1e9, 1e-9*f_full(length, young['full'][0]), '-r')
plt.legend(["Datos","Ajuste"])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('Full_Fit', name), bbox_inches='tight')

#%% 2G) FREQUENCY AND LENGTH
# --> Try a nonlinear fit using the Andrea model directly

young['andrea'] = {}
factor['andrea'] = {}
chi_squared['andrea'] = {}

initial_guess = (64e9, .2)
rsq, parameters = iva.nonLinearFit(length, frequency, f_andrea, 
                                   initial_guess=initial_guess, 
                                   bounds=([1e9,0], [np.infty, 1]), 
                                   showplot=False)
young['andrea'] = parameters[0]
factor['andrea'] = parameters[1]
del parameters
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young['andrea'][0], 
                                                        young['andrea'][1], 
                                                        units="Pa")))
print(r"Factor porcentual: {}%".format(ivu.errorValueLatex(
        factor['andrea'][0]*100,
        factor['andrea'][1]*100)))

chi_squared['andrea'] = sum( (f_andrea(length, 
                                       young['andrea'][0], 
                                       factor['andrea'][0]) 
                            - frequency)**2 ) 
chi_squared['andrea'] = chi_squared['andrea'] / len(length)

plt.figure()
ax = plt.subplot()
plt.title('Ajuste completo aproximado con factor')
plt.plot(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.plot(length*1e9, 1e-9*f_andrea(length, young['andrea'][0], 
                                     factor['andrea'][0]), '-r')
plt.legend(["Datos","Ajuste"])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('Andrea_Fit', name), bbox_inches='tight')

#%% 2H) FREQUENCY AND LENGTH
# --> Try a nonlinear fit directly using the complex model with free factor

young['iv'] = {}
factor['iv'] = {}
chi_squared['iv'] = {}
    
initial_guess = (64e9, .2)
rsq, parameters = iva.nonLinearFit(length, frequency, f_iv, 
                                   initial_guess=initial_guess, 
                                   bounds=([1e9,0], [np.infty, 1]), 
                                   showplot=False)
young['iv'] = parameters[0]
factor['iv'] = parameters[1]
del parameters
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young['iv'][0], 
                                                        young['iv'][1], 
                                                        units="Pa")))
print(r"Factor porcentual: {}%".format(ivu.errorValueLatex(
        factor['iv'][0]*100,
        factor['iv'][1]*100)))

chi_squared['iv'] = sum( (f_iv(length, young['iv'][0], factor['iv'][0]) 
                          - frequency)**2 ) 
chi_squared['iv'] = chi_squared['iv'] / len(length)
  
plt.figure()
ax = plt.subplot()
plt.title('Ajuste completo con factor')
plt.plot(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.plot(length*1e9, 1e-9*f_andrea(length, young['iv'][0], 
                                     factor['iv'][0]), '-r')
plt.legend(["Datos","Ajuste"])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('IV_Fit', name), bbox_inches='tight')

#%% 2I) FREQUENCY AND LENGTH
# --> Try a nonlinear fit using the complex model with free K1, K2, Young

#rsq, parameters = iva.nonLinearFit(length, frequency, f_free, 
#                                   initial_guess=(Young, K1, K2),
#                                   showplot=False)
#young['free'] = parameters[0]
#k1['free'] = parameters[1]
#k2['free'] = parameters[2]
#del parameters
#print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young['free'][0], 
#                                                        young['free'][1], 
#                                                        units="Pa")))
#
#chi_squared['free'] = sum( (f_free(length, young['free'][0], k1['free'][0], k2['free'][0]) 
#                               - frequency)**2 ) 
#chi_squared['free'] = chi_squared['free'] / len(length)
#
#plt.figure()
#ax = plt.subplot()
#plt.title('Ajuste modelo completo libre')
#plt.plot(length*1e9, frequency*1e-9,'o')
#plt.ylabel('Frecuencia (GHz)')
#plt.xlabel('Longitud (nm)')
#plt.plot(length*1e9, 1e-9*f_free(length, young['free'][0], k1['free'][0], k2['free'][0]), '-r')
#plt.legend(["Datos","Ajuste completo"])
#ax.minorticks_on()
#ax.tick_params(axis='y')
#ax.tick_params(axis='y', which='minor', length=0)
#ax.grid(axis='both', which='both')
#plt.show()
#
## Save plot
#plt.savefig(figsFilename('Free_Fit', name), bbox_inches='tight')

#%% 2*) FREQUENCY AND LENGTH
# --> Transform values into a more fast-forward way of reading them

young_str = {}
factor_str = {}
chi_squared_str = {}

for typ in ['simple', 'full', 'andrea', 'iv']:
    young_str[typ] = ivu.errorValueLatex(young[typ][0], 
                                     young[typ][1], 
                                     units="Pa")
    try:
        factor_str[typ] = ivu.errorValueLatex(factor[typ][0]*100, 
                                             factor[typ][1]*100)
    except:
        typ
    chi_squared_str[typ] = '{:.2E}'.format(chi_squared[typ])

#%% 2**) FREQUENCY AND LENGTH
# --> Final plots of fits

plt.figure()
ax = plt.subplot()
plt.loglog(length*1e9, frequency*1e-9, 'o', label='Datos')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.loglog(length*1e9, 1e-9*f_simple(length, young['simple'][0]), '-k', 
           label='Ajuste modelo simple')
plt.loglog(length*1e9, 1e-9*f_full(length, young['full'][0]), '-r', 
           label='Ajuste modelo completo')
plt.loglog(length*1e9, 1e-9*f_iv(length, young['iv'][0], factor['iv'][0]), '-g', 
           label=r'Ajuste modelo completo con factor')
plt.loglog(length*1e9, 1e-9*f_andrea(length, young['andrea'][0],
                                     factor['andrea'][0]), '--c', 
           label=r'Ajuste modelo completo aproximado con factor')
plt.legend()
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('LoglogFinal', name), bbox_inches='tight')


plt.figure()
ax = plt.subplot()
plt.plot(length*1e9, frequency*1e-9, 'o', label='Datos')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.plot(length*1e9, 1e-9*f_simple(length, young['simple'][0]), '-k', 
           label='Ajuste modelo simple')
plt.plot(length*1e9, 1e-9*f_full(length, young['full'][0]), '-r', 
           label='Ajuste modelo completo')
plt.plot(length*1e9, 1e-9*f_iv(length, young['iv'][0], factor['iv'][0]), '-g', 
           label=r'Ajuste modelo completo con factor')
plt.plot(length*1e9, 1e-9*f_andrea(length, young['andrea'][0],
                                     factor['andrea'][0]), '--c', 
           label=r'Ajuste modelo completo aproximado con factor')
plt.legend()
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('FitsFinal', name), bbox_inches='tight')

#%% 2***) FREQUENCY AND LENGTH
# --> Final plots of F vs L

# Data for predictions
young_predict_select = 64e9
factor_predict = [0, .1, .2] # fraction that represents bound

# Theory predictions
frequency_predict = np.array([f_iv(length, young_predict_select, c) 
                              for c in factor_predict]).T

# Make a plot with fit and predictions
plt.figure()
ax = plt.subplot()
plt.title('Análisis de resultados')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.plot(length*1e9, frequency*1e-9,'o')
plt.plot(length*1e9, 1e-9*f_simple(length, young['simple'][0]), '-r')
for freq in frequency_predict.T: plt.plot(length*1e9, freq*1e-9, ':')
del freq
plt.legend(["Datos", r"Ajuste en vacío con {}".format(young_str['simple'])] + 
          ["Predicción inmerso con {:.0f} GPa al {:.0f}%".format(
                  young_predict_select/1e9,
                  f*100) for f in factor_predict])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('PlotSuperFinal', name), bbox_inches='tight')

# Make a loglog plot with fit and predictions
plt.figure()
ax = plt.subplot()
plt.title('Análisis de resultados')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.loglog(length*1e9, frequency*1e-9,'o')
plt.loglog(length*1e9, 1e-9*f_simple(length, young['simple'][0]), '-r')
for freq in frequency_predict.T: plt.loglog(length*1e9, freq*1e-9, ':')
del freq
plt.legend(["Datos", r"Ajuste en vacío con {}".format(young_str['simple'])] + 
          ["Predicción inmerso con {:.0f} GPa al {:.0f}%".format(
                  young_predict_select/1e9,
                  f*100) for f in factor_predict])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()

# Save plot
plt.savefig(figsFilename('LoglogSuperFinal', name), bbox_inches='tight')

del factor_predict, young_predict_select, frequency_predict

"""
Multiple legends on the same Axes
https://matplotlib.org/users/legend_guide.html
"""

#%% 3) FREQUENCY VS Q AND WIDTH

# Plot results 
fig, ax1 = plt.subplots()

# Frequency vs width plot, lower axis
ax1.set_xlabel('Ancho (nm)', color='tab:red')
ax1.set_ylabel('Frecuencia (GHz)')
ax1.plot(sem_data[:,0], fits_data[:,0], 'ro')
ax1.tick_params(axis='x', labelcolor='tab:red')

# Frequency vs quality factor, upper axis
ax2 = ax1.twiny()  # Second axes that shares the same y-axis
ax2.set_xlabel('Factor de calidad (u.a.)', color='tab:blue')
ax2.plot(sem_data[:,4], fits_data[:,0], 'bx')
ax2.tick_params(axis='x', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
ax1.grid(axis='both')

# Save plot
plt.savefig(figsFilename('FvsWyQ', name), bbox_inches='tight')

"""LIGO1: La frec no parece verse influenciada por ninguna de las dos variables."""

#%% 4) DAMPING TIME VS LENGTH AND WIDTH

# Plot results 
fig, ax1 = plt.subplots()

# Frequency vs width plot, lower axis
ax1.set_xlabel('Ancho (nm)', color='tab:red')
ax1.set_ylabel('Tiempo de decaimiento (GHz)')
ax1.plot(sem_data[:,0], fits_data[:,2], 'ro')
ax1.tick_params(axis='x', labelcolor='tab:red')

# Frequency vs quality factor, upper axis
ax2 = ax1.twiny()  # Second axes that shares the same y-axis
ax2.set_xlabel('Longitud (nm)', color='tab:blue')
ax2.plot(sem_data[:,2], fits_data[:,2], 'bx')
ax2.tick_params(axis='x', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
ax1.grid(axis='both')

# Save plot
plt.savefig(figsFilename('TauvsLyW', name), bbox_inches='tight')

"""LIGO1: El tiempo de decaimiento a lo sumo decae con ambas. Pero el modelo 
predice que aumenta con L**2 Oo"""

#%% 5) DAMPING TIME VS LENGTH

# Fit
viscosity = iva.nonLinearFit(length, damping_time, tau_simple, 
                             showplot=False)[1][0]

# Plot
plt.figure()
ax = plt.subplot()
plt.plot(length, damping_time, '.')
plt.plot(length, tau_simple(length, viscosity[0]), '-r')
plt.xlabel('Longitud $L$ (m)')
plt.ylabel(r'Tiempo de decaimiento $\tau$ (s)')
plt.title(r'Ajuste')
plt.legend(['Datos', 'Ajuste'])
plt.grid(axis='x')

print(r"Viscosidad: {}".format(ivu.errorValueLatex(
    viscosity[0], 
    viscosity[1], 
    units="Pa.s")))

#%% *) HISTOGRAMS

fig = plt.figure()
grid = plt.GridSpec(1, 2, wspace=0)

ax = plt.subplot(grid[0,0])
bins_limits = ax.hist(sem_data[:,2])[1]
plt.xlabel("Longitud (nm)")
plt.ylabel("Repeticiones")

#mean_frequencies = []
#for Fi, Ff in zip(bins_limits[:-1], bins_limits[1:]):
#    mean_frequencies = np.mean(fits_data[Fi<=fits_data[:,0]<Ff,0])

ax = plt.subplot(grid[0,1])
bins_limits = ax.hist(fits_data[:,0])[1]
plt.xlabel("Frecuencia (GHz)")

# Save plot
plt.savefig(figsFilename('Hists', name), bbox_inches='tight')

#mean_frequencies = []
#for Fi, Ff in zip(bins_limits[:-1], bins_limits[1:]):
#    mean_frequencies = np.mean(fits_data[Fi<=fits_data[:,0]<Ff,0])

#%% *) BOX PLOTS

fig = plt.figure()
grid = plt.GridSpec(2, 2, wspace=0.5, hspace=0)

ax = plt.subplot(grid[0,0])
ax.boxplot(fits_data[:,0], showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2})
for l in ax.get_xticklabels():
    l.set_visible(False)
del l
plt.ylabel("Frecuencia (GHz)")
ax.tick_params(axis='y', direction='in')

ax = plt.subplot(grid[0,1])
ax.boxplot(sem_data[:,2], showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2})
for l in ax.get_xticklabels():
    l.set_visible(False)
del l
plt.ylabel("Longitud (nm)")
ax.tick_params(axis='y', direction='in')

ax = plt.subplot(grid[1,1])
ax.boxplot(sem_data[:,0], showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2})
for l in ax.get_xticklabels():
    l.set_visible(False)
del l
plt.ylabel("Ancho (nm)")
ax.tick_params(axis='y', direction='in')

ax = plt.subplot(grid[1,0])
ax.boxplot(fits_data[:,2], showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2})
plt.ylabel("Factor de calidad")
ax.tick_params(axis='y', direction='in')

# Save plot
plt.savefig(figsFilename('Boxplots', name), bbox_inches='tight')


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

#%% *) EXTRA CODE
    
ivu.copy(ivu.errorValueLatex(np.mean(sem_data[:,0]),
                             max(np.mean(sem_data[:,1]), np.std(sem_data[:,0])),
                             units='nm'))

ivu.copy(ivu.errorValueLatex(np.mean(fits_data[:,2]),
                             np.std(fits_data[:,2]),
                             units='GHz'))

ivu.copy(ivu.errorValueLatex(np.mean(params[:,1]),
                             np.std(params[:,1]),
                             units='$\mu$W'))