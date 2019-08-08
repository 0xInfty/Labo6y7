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
"""
rods_filename = os.path.join(home, r'Análisis\Rods_LIGO1.txt')
sem_series = ['LIGO1_1']
sem_short_series = lambda series : '{}'#series.split('_')[1]+' {}'
name = 'LIGO1'
"""

rods_filename = os.path.join(home, r'Análisis\Rods_M135.txt')
sem_series = ['M135_5_1D', 'M135_7B_1D']
sem_short_series = lambda series : series.split('_')[1]+' {}'
name = 'M135'

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
diameter = 27e-9 # m for rods
midlength = 85e-9 # m for rods
viscosity = 2e-3 # Pa.s for gold
Young = np.mean([71.2e9, 74.8e9])  # Pa for fused silica
density_s = np.mean([2.17e3, 2.22e3]) # kg/m3 for fused silica
K1 = Shear * 2.75 # Pa
K2 = np.pi * diameter * np.sqrt(density_s * Shear)

# Space to save results
young = {}
factor = {}
chi_squared = {}
area_mode = ['circ', 'rect']

# Theory models
def f_simple(length, young):
    f_0 = (np.sqrt(young/density) / (2 * length))
    return f_0

def f_mid(length, young):
    f_0 = f_simple(length, young)
    beta = ( viscosity / (length * density) )**2 / 2
    f = np.sqrt(f_0**2 - beta**2 / 4)
    return f

def f_area_definer(mode):

    if mode=='rect':
        area = np.pi * diameter * midlength
    elif mode=='circ':
        area = np.pi * (diameter**2) / 4
    else:
        raise ValueError("Wrong key!")

    def f_andrea(length, young, factor=1):
        f_0 = f_simple(length, young)
        f = np.sqrt( (2*np.pi*f_0)**2 + factor*K1/(density*area) ) 
        f = f /(2*np.pi)
        return f

    def f_complex(length, young):
        f_0 = f_simple(length, young)
        beta = ( viscosity / (length * density) )**2 / 2
        K1_term = K1 / ( np.pi**2 * density * area )
        K2_subterm = K2 / ( 2 * np.pi * density * area )
        f = np.sqrt(f_0**2 + K1_term/4 - (K2_subterm + beta)**2/4 )
        return f
    
    def f_vall(length, young, factor=1):
        f_0 = f_simple(length, young)
        beta = ( viscosity / (length * density) )**2 / 2
        K1_term = K1 / ( np.pi**2 * density * factor * area )
        K2_subterm = K2 / ( 2 * np.pi * density * factor * area )
        f = np.sqrt(f_0**2 + K1_term/4 - (K2_subterm + beta)**2/4 )
        return f

    return f_andrea, f_complex, f_vall

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
damping_time = fits_data[:,1]*1e12 # s
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

# Save plot
plt.savefig(figsFilename('FyQvsRod', name), bbox_inches='tight')

#%% 2A) FREQUENCY AND LENGTH
# --> Try out some known values

f_andrea, f_complex, f_vall = f_area_definer('circ')

# Data
young = [45e9, 64e9, 78e9]
factor = [0, .1, .2, .3, .4] # fraction that represents bound

# Theory predictions
f_0 = np.array([f_simple(l, y) for l in length for y in young])
f_0 = f_0.reshape([len(length), len(young)])
f = np.array([f_andrea(l, young[-1], c) for l in length for c in factor])
f = f.reshape([len(length), len(factor)])

# Make a plot for the simpler model
plt.figure()
ax = plt.subplot()
plt.title('Modelo simple')
plt.loglog(length*1e9, frequency*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in f_0.T: plt.loglog(length*1e9, freq*1e-9, '-')
del freq
plt.legend(["Datos"] + ["{} GPa".format(y/1e9) for y in young])
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
plt.show()
for l in ax.get_xticklabels():
    l.set_visible(False)
del l

# Make one too for the complex model
plt.figure()
plt.title('Modelo complejo con Young {} GPa'.format(young[-1]/1e9))
plt.loglog(length*1e9, fits_data[:,0]*1e-9,'o')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
for freq in f.T: plt.loglog(length*1e9, freq*1e-9, '-')
del freq
plt.legend(["Datos"] + ["Sumergido {:.0f}%".format(c*100) for c in factor])

"""LIGO1: Decidimos que esto no es necesario si hacemos ajustes"""

#%% 2B) FREQUENCY AND LENGTH
# --> Try a linear fit on f vs 1/L, which corresponds to the simple model

#rsq, m, b = iva.linearFit(1/length, fits_data[:,0], M = True)
#plt.ylabel('Frecuencia (GHz)')
#plt.xlabel('Inverso de longitud (1/nm)')
#young_fit = 4 * density * (m[0]**2)
#young_fit_error = np.abs(8 * density * m[1] * m[0])
#print(r"Módulo de Young: {}".format(ivu.errorValueLatex(young_fit, 
#                                                        young_fit_error, 
#                                                        units="Pa")))

"""LIGO 1: Da pendiente -1.5 aproximadamente, un sinsentido físico"""

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

"""LIGO1
Módulo de Young: (78.0$\pm$3.2) GPa
Chi Squared: 1.0990569627612585e+18
Hermoso. Absolutamente hermoso :)
"""

#%% 2E) FREQUENCY AND LENGTH
# --> Try a nonlinear fit directly using the mid model

young['mid'] = iva.nonLinearFit(length, frequency, 
                                f_mid, showplot=False)[-1][0]
print(r"Módulo de Young: {}".format(ivu.errorValueLatex(
            young['mid'][0], 
            young['mid'][1], 
            units="Pa")))
chi_squared['mid'] = sum( (f_mid(length, young['mid'][0]) - frequency)**2 ) 
chi_squared['mid'] = chi_squared['mid'] / len(length)

"""LIGO1
Módulo de Young: (78.0$\pm$3.2) GPa
Chi Squared 1.0517810075757833e+18
Como la viscosidad del oro es muy pequeña y las escalas de longitud 
son nanométricas, el término beta se puede despreciar frente a la frecuencia 
omega_0. El máximo al que llega beta**2/4 es 20 órdenes menor al de f_0**2.
"""
#%% 2F) FREQUENCY AND LENGTH
# --> Try a nonlinear fit using the Andrea model directly

young['andrea'] = {}
factor['andrea'] = {}
chi_squared['andrea'] = {}

for a in area_mode:

    f_andrea, f_complex, f_vall = f_area_definer(a)
    
    initial_guess = (Young, .2)
    rsq, parameters = iva.nonLinearFit(length, frequency, f_andrea, 
                                       initial_guess=initial_guess, 
                                       bounds=([1e9,0], [np.infty, 1]), 
                                       showplot=False)
    y = parameters[0]
    f = parameters[1]
    print(r"Módulo de Young: {}".format(ivu.errorValueLatex(y[0], 
                                                            y[1], 
                                                            units="Pa")))
    print(r"Factor porcentual: {}%".format(ivu.errorValueLatex(f[0]*100,
                                                               f[1]*100)))
    
    ch = sum( (f_andrea(length, y[0], f[0]) - frequency)**2 ) 
    ch = ch / len(length)

    young['andrea'][a] = y
    factor['andrea'][a] = f
    chi_squared['andrea'][a] = ch

del a, y, f, ch

"""LIGO1

Si las CI son (64e9, .2)... Usando área transversal circular...
Módulo de Young: (78.4$\pm$24.1) GPa
Factor porcentual: (6.4$\pm$216000000000000.0)$10^{-13}$%
Chi Squared: 1.051781007575788e+18

Usando área transversal rectangular...
Módulo de Young: (78.0$\pm$24.4) GPa
Factor porcentual: (2.0$\pm$2770000000000000.0)$10^{-13}$%
Chi Squared: 1.0517810075757824e+18

Todo indica que el algoritmo quiere bajar lo más que puede el factor.
"""

#%% 2G) FREQUENCY AND LENGTH
# --> Try a nonlinear fit directly using the complex model

young['complex'] = {}
chi_squared['complex'] = {}

for a in area_mode:

    f_andrea, f_complex, f_vall = f_area_definer(a)

    
    initial_guess = (Young)
    rsq, y = iva.nonLinearFit(length, frequency, f_complex, 
                              initial_guess=initial_guess, 
                              bounds=([1e9], [np.infty]), 
                              showplot=False)
    y = parameters[0]
    print(r"Módulo de Young: {}".format(ivu.errorValueLatex(y[0], 
                                                            y[1], 
                                                            units="Pa")))
    
    ch = sum( (f_complex(length, y[0]) - frequency)**2 ) 
    ch = ch / len(length)

    young['complex'][a] = y
    chi_squared['complex'][a] = ch

del a, y, ch

"""LIGO1

Usando área transversal circular...
Módulo de Young: (78.4$\pm$24.1) GPa
Chi Squared: 3.5884968302104535e+19

Usando área transversal rectangular...
Módulo de Young: (78.4$\pm$24.1) GPa
Chi Squared: 1.527528834987675e+18
"""

#%% 2H) FREQUENCY AND LENGTH
# --> Try a nonlinear fit directly using the complex model with free factor

young['vall'] = {}
factor['vall'] = {}
chi_squared['vall'] = {}

for a in area_mode:

    f_andrea, f_complex, f_vall = f_area_definer(a)
    
    initial_guess = (Young, .2)
    rsq, parameters = iva.nonLinearFit(length, frequency, f_vall, 
                                       initial_guess=initial_guess, 
                                       bounds=([1e9,0], [np.infty, 1]), 
                                       showplot=False)
    y = parameters[0]
    f = parameters[1]
    print(r"Módulo de Young: {}".format(ivu.errorValueLatex(y[0], 
                                                            y[1], 
                                                            units="Pa")))
    print(r"Factor porcentual: {}%".format(ivu.errorValueLatex(f[0]*100,
                                                               f[1]*100)))
    
    ch = sum( (f_vall(length, y[0], f[0]) - frequency)**2 ) 
    ch = ch / len(length)

    young['vall'][a] = y
    factor['vall'][a] = f
    chi_squared['vall'][a] = ch
  
del a, y, f, ch
    
"""LIGO1

Si las CI son (Young, 1)... Usando el área transversal circular...
Módulo de Young: (73.0$\pm$139.0) GPa
Factor porcentual: (14.2$\pm$3.1)%
Chi Squared: 1.8031943401351125e+19

Usando área transversal rectangular...
Módulo de Young: (69.8$\pm$25.1) GPa
Factor porcentual: (100.0$\pm$291.0)%
Chi Squared: 1.1227733049602669e+18
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
for l in ax.get_yticklabels():
    l.set_visible(False)
del l

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
for l in ax.get_xticklabels():
    l.set_visible(False)
del l
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

for typ in ['andrea','complex','vall']:
    for a in ['circ', 'rect']:
        young[typ][a] = ivu.errorValueLatex(young[typ][a][0], 
                                            young[typ][a][1], 
                                            units="Pa")
        try:
            factor[typ][a] = ivu.errorValueLatex(factor[typ][a][0]*100, 
                                                 factor[typ][a][1]*100)
        except:
            a
        chi_squared[typ][a] = '{:.2e}'.format(chi_squared[typ][a])
for typ in ['simple', 'mid']:
    young[typ] = ivu.errorValueLatex(young[typ][0], 
                                     young[typ][1], 
                                     units="Pa")
    chi_squared[typ] = '{:.2E}'.format(chi_squared[typ])

#%% *) EXTRA CODE
    
ivu.copy(ivu.errorValueLatex(np.mean(sem_data[:,0]),
                             max(np.mean(sem_data[:,1]), np.std(sem_data[:,1])),
                             units='nm'))

ivu.copy(ivu.errorValueLatex(np.mean(fits_data[:,2]),
                             np.std(fits_data[:,2]),
                             units='GHz'))

ivu.copy(ivu.errorValueLatex(np.mean(params[:,3]),
                             np.std(params[:,3]),
                             units='nm'))