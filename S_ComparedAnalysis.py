# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:58:43 2019

@author: Vall
"""

import numpy as np
#import matplotlib.pyplot as plt
import os
import iv_save_module as ivs
#import iv_utilities_module as ivu
#import iv_analysis_module as iva

#%% PARAMETERS ----------------------------------------------------------------

# Main folder's path
home = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7'

# For each data to compare, we need one value on each list
desired_frequency = [12, 17] # in GHz
series = ['LIGO1_PostUSA', 'LIGO1']
sem_series = ['LIGO1_1', 'LIGO1_1']

# Some functions and variables to manege filenames
def semFilename(sem_series, home=home):
    """Given a series 'M135_7B_1D', returns path to SEM data"""
    filename = 'Resultados_SEM_{}.txt'.format(series)
    sem_series = sem_series.split('_') # From 'M_20190610_01' take '20190610'
    sem_filename = os.path.join(home, 'Muestras\SEM', *sem_series, filename)
    return sem_filename
rodsFilename = lambda series : os.path.join(home, 
                                           r'Análisis\Rods_{}.txt'.format(
                                                   series))
paramsFilename = lambda series : os.path.join(home, 
                                              r'Análisis/Params_{}.txt'.format(
                                                   series))
figsFilename = lambda fig_name : os.path.join(home, fig_name+'.png')
figs_folder = r'Análisis/ComparedAnalysis'+series
figs_extension = '.png'

#%% LOAD DATA -----------------------------------------------------------------

filenames = []
rods = []
params = []
fits_data = []
fits_footer = []

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
    index = [params_filenames.index(f) for f in filenames]
    sparams = sparams[index,:]
    params_header = ['Amplitud (mVpp)', 'Potencia Pump post-MOA (muW)', 
                     'Longitud de onda (nm)', 'Ancho medio de la campana (nm)']
    del params_filenames, index, amplitude, power, wavelength, spectral_width
    
    # Now create a list of folders for each filename    
    fits_filenames = [ivs.filenameToFitsFilename(file, home) for file in filenames]
    
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
    for rod, fit in zip(rods, sfits_data):
        try:
            i = np.argmin(abs(fit[:,0] - desired_frequency*np.ones(fit.shape[0])))
            fits_new_data.append([*fit[i,:]])
        except IndexError:
            fits_new_data.append([*fit])
    sfits_data = np.array(fits_new_data)
    frequency = sfits_data[:,0]*1e9 # Hz
    damping_time = sfits_data[:,1]*1e-12 # s
    quality_factor = sfits_data[:,2]
    del rod, fit, i, fits_new_data
    
    #### ACÁ ME QUEDEÉEEEEEEEE
    
    # Also create a list of folders for SEM filenames
    sem_filenames = [filenameToSEMFilename(s) for s in sem_series]
    
    # Then load data from SEM dimension analysis
    sem_data = []
    sem_footer = []
    sem_rods = []
    for f in zip(sem_filenames, sem_series):
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