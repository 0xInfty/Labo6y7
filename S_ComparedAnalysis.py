# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:58:43 2019

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
load_sem = True
fit_series = 0

# For each data to compare, we need one value on each list
desired_frequency = [12, 17] # in GHz
full_series = ['Aire', 'Ta2O5']
series = ['LIGO1', 'LIGO1_PostUSA']
sem_series = ['LIGO1_1', 'LIGO1_1']
name = 'FusedSilica'

# Some functions and variables to manege filenames
def semFilename(sem_series, home=home):
    """Given a series 'M135_7B_1D', returns path to SEM data"""
    filename = 'Resultados_SEM_{}.txt'.format(sem_series)
    sem_series = sem_series.split('_') # From 'M_20190610_01' take '20190610'
    sem_filename = os.path.join(home, 'Muestras\SEM', *sem_series, filename)
    return sem_filename
rodsFilename = lambda series : os.path.join(home, 
                                           r'Análisis\Rods_{}.txt'.format(
                                                   series))
paramsFilename = lambda series : os.path.join(home, 
                                              r'Análisis\Params_{}.txt'.format(
                                                   series))
figsFilename = lambda fig_name : os.path.join(home, fig_name+'.png')
figs_folder = r'Análisis/ComparedAnalysis_{}'.format(name)
figs_extension = '.png'
overwrite = True

#%% MODELS --------------------------------------------------------------------

# Physics' Parameters
physics = {}
physics['density_gold'] = 19.3e3 # kg/m3 for gold
physics['shear_fsilica'] = np.mean([30.8e9, 32.3e9]) # Pa for fused silica
physics['diameter'] = 27.7e-9 # m for rods
physics['midlength'] = 85e-9 # m for rods
physics['viscosity_gold'] = 2e-3 # Pa/s for gold
physics['young_gold'] = np.mean([71.2e9, 74.8e9])  # Pa for fused silica
physics['density_fsilica'] = np.mean([2.17e3, 2.22e3]) # kg/m3 for fused silica
physics['area'] = np.pi * physics['diameter']**2 / 4
physics['K1'] = physics['shear_fsilica'] * 2.75 # Pa
physics['K2'] = np.pi * physics['diameter'] 
aux = np.sqrt(physics['density_fsilica'] * physics['shear_fsilica']) 
physics['K2'] = physics['K2'] * aux # Pa.s (viscosity's units)
del aux

# Theory models
def f_simple(length, young):
    f_0 = (np.sqrt(young/physics['density_gold']) / (2 * length))
    return f_0

#%% LOAD DATA -----------------------------------------------------------------

filenames = []
rods = []
params = []
fits_data = []
fits_footer = []
frequency = []
damping_time = []
quality_factor = []
if load_sem:
    sem_data = []
    sem_footer = []
    length = []
    width = []

for s, ss, f in zip(series, sem_series, desired_frequency):

    # Look for the list of rods and filenames
    sfilenames = [] # Will contain filenames like 'M_20190610_01'
    srods = [] # Will contain rods' positions like '1,2'
    with open(rodsFilename(s), 'r') as file:
        for line in file:
            if line[0]!='#':
                sfilenames.append(line.split('\t')[0]) # Save filenames
                aux = line.split('\t')[1:]
                aux = r' '.join(aux)
                srods.append(aux.split('\n')[0]) # Save rods
                del aux
        del line
    
    # Then load parameters
    params_filenames = [] # Will contain filenames like 'M_20190610_01'
    amplitude = []
    power = []
    wavelength = []
    spectral_width = []
    with open(paramsFilename(s), 'r') as file:
        for line in file:
            if line[0]!='#':
                params_filenames.append(line.split('\t')[0])
                amplitude.append(float(line.split('\t')[1]))
                power.append(float(line.split('\t')[2]))
                wavelength.append(float(line.split('\t')[3]))
                spectral_width.append(float(line.split('\t')[-1]))
        del line
    sparams = np.array([amplitude, power, wavelength, spectral_width]).T
    index = [params_filenames.index(f) for f in sfilenames]
    sparams = sparams[index,:]
    params_header = ['Amplitud (mVpp)', 'Potencia Pump post-MOA (muW)', 
                     'Longitud de onda (nm)', 'Ancho medio de la campana (nm)']
    del params_filenames, index, amplitude, power, wavelength, spectral_width
    
    # Now create a list of folders for each filename    
    fits_filenames = [ivs.filenameToFitsFilename(file, home) for file in sfilenames]
    
    # Load data from each fit
    sfits_data = []
    sfits_footer = []
    for file in fits_filenames:
        data, fits_header, footer = ivs.loadTxt(file)
        sfits_data.append(data)
        sfits_footer.append(footer)
    del file, data, footer, fits_filenames
    
    # Keep only the fit term that has the closest frequency to the desired one
    fits_new_data = []
    for rod, fit in zip(srods, sfits_data):
        try:
            i = np.argmin(abs(fit[:,0] - f*np.ones(fit.shape[0])))
            fits_new_data.append([*fit[i,:]])
        except IndexError:
            fits_new_data.append([*fit])
    sfits_data = np.array(fits_new_data)
    sfrequency = sfits_data[:,0]*1e9 # Hz
    sdamping_time = sfits_data[:,1]*1e-12 # s
    squality_factor = sfits_data[:,2]
    del rod, fit, i, fits_new_data
    
    if load_sem:
    
        # Also load data from SEM dimension analysis
        ssem_data, sem_header, ssem_footer = ivs.loadTxt(semFilename(ss))
    #    other_data = []
    #    for d in ssem_data:
    #        for di in d:
    #            other_data.append([*di])
    #    other_data =  np.array(other_data)
    #    del d, di
        
        # Now lets put every rod in the same order for SEM and fits
        index = [ssem_footer['rods'].index(r) for r in srods]
        ssem_data = np.array([ssem_data[i,:] for i in index])
        slength = ssem_data[:,2] * 1e-9 # m
        swidth = ssem_data[:,0] * 1e-9 # m
        del index
    
#    # Now we can filter the results
#    index = np.argsort(sfrequency) # Remove the two lowest frequencies
#    if load_sem:
#       slength = slength[index[2:]]
#       swidth = swidth[index[2:]]
#    sfrequency = sfrequency[index[2:]]
#    sdamping_time = sdamping_time[index[2:]]
#    squality_factor = squality_factor[index[2:]]
#    del index
    
    # Since I'll be analysing frequency vs length mostly...
    if load_sem:
        index = np.argsort(slength)
        slength = slength[index]
        swidth = swidth[index]
        sfrequency = sfrequency[index]
        sdamping_time = sdamping_time[index]
        squality_factor = squality_factor[index]
        del index
    
    # Now add all that data to a list outside the loop
    filenames.append(sfilenames)
    rods.append(srods)
    params.append(sparams)
    fits_data.append(sfits_data)
    fits_footer.append(sfits_footer)
    frequency.append(sfrequency)
    damping_time.append(sdamping_time)
    quality_factor.append(squality_factor)
    if load_sem:
        sem_data.append(ssem_data)
        sem_footer.append(ssem_footer)
        length.append(slength)
        width.append(swidth)
    
del sfilenames, srods, sfits_data, sfits_footer, sfrequency, sdamping_time
del squality_factor, ssem_data, ssem_footer, slength, swidth

#%%

"""

### HASTA ACÁ LLEGUÉ

tables = []

for i, s, fs in zip(range(len(series)), series, full_series):
    
    # Prepare important data for a table 
    items = []
    for i in range(len(rods)):
        h = '\t'.join(ivu.errorValue(sem_dataDSAGASDGGD[i,2], sem_data[i,3]))
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
    tables.append('\n'.join([heading, items]))
    del heading, items
    
    # Save all important data to a single file
    whole_filename = os.path.join(home, r'Análisis/Resultados_Totales_{}.txt'.format(series))
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
   
"""
#%% *) FREQUENCY PER ROD

# Plot results for the different rods
fig, ax1 = plt.subplots()

# Frequency plot, right axis
ax1.set_xlabel('Antena')
ax1.set_ylabel('Frecuencia (GHz)', color='tab:red')
ax1.plot(fits_data[0][:,0], 'ro', fits_data[1][:,0], 'rx')
ax1.tick_params(axis='y', labelcolor='tab:red')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
ax1.legend(full_series)
plt.show()

# Format graph
plt.xticks(np.arange(len(rods[0])), rods[0], rotation='vertical')
plt.grid(which='both', axis='x')
ax1.tick_params(length=2)
ax1.grid(axis='x', which='minor')
ax1.tick_params(axis='x', labelrotation=90)
plt.show()

# Save plot
ivs.saveFig(figsFilename('FvsRod'), extension=figs_extension, 
            folder=figs_folder, overwrite=overwrite)

#%% *) BOXPLOTS

fig = plt.figure()
grid = plt.GridSpec(1, 2, wspace=0.5, hspace=0)

ax = plt.subplot(grid[0,0])
ax.boxplot([fits_data[0][:,0], fits_data[1][:,0]], 
           showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2},
           labels=full_series)
plt.ylabel("Frecuencia (GHz)")
ax.tick_params(axis='y', direction='in')
    
ax = plt.subplot(grid[0,1])
ax.boxplot([fits_data[0][:,2], fits_data[1][:,2]], 
           showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2},
           labels=full_series)
plt.ylabel("Factor de calidad")
ax.tick_params(axis='y', direction='in')

ivs.saveFig(figsFilename('Boxplots'), extension=figs_extension, 
            folder=figs_folder, overwrite=overwrite)

#%% *) FREQUENCY AND LENGTH FIT

if load_sem:
    young_gold = iva.nonLinearFit(length[fit_series], frequency[fit_series], 
                                  f_simple, showplot=False)[-1][0]
    print(r"Módulo de Young: {}".format(ivu.errorValueLatex(*young_gold, 
                                                            units="Pa")))
    chi_squared = sum( (f_simple(length[fit_series],
                                 young_gold[0]) - frequency[fit_series])**2 ) 
    chi_squared = chi_squared / len(length[fit_series])

#%% *) FREQUENCY AND LENGTH
# --> Final plots of F vs L

# Make a plot with fit and predictions
plt.figure()
ax = plt.subplot()
plt.title('Análisis de resultados')
plt.ylabel('Frecuencia (GHz)')
plt.xlabel('Longitud (nm)')
plt.plot(length[0]*1e9, frequency[0]*1e-9,'ok', label=full_series[0])
plt.plot(length[1]*1e9, frequency[1]*1e-9, 'xb', label=full_series[1])
plt.plot(length[fit_series]*1e9, 
         1e-9*f_simple(length[fit_series], young_gold[0]), 
         '-r', 
         label=r"Ajuste en vacío con {}".format(ivu.errorValueLatex(
                 *young_gold, 
                 units="Pa")))
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='both', which='both')
ax.legend()
plt.show()

# Save plot
ivs.saveFig(figsFilename('PlotSuperFinal'), extension=figs_extension, 
            folder=figs_folder, overwrite=overwrite)

"""
Multiple legends on the same Axes
https://matplotlib.org/users/legend_guide.html
"""