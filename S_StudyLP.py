# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:19:29 2019

@author: Vall
"""

import iv_analysis_module as iva
import matplotlib.pyplot as plt
import iv_save_module as ivs
import iv_utilities_module as ivu
import numpy as np
import os

#%%

# Parameters
home = r'C:\Users\Luciana\OneDrive\Labo 6 y 7'
names = ['M_20190605_11']
desired_frequency = 9 # Desired frequency for the ideal fit
Ni = 40 # How many index around the main one we'll try for the initial time

#%%

# Data to collect while iterating
jmean = [] # Mean index
jgood = [] # Index that allow fitting
jreallygood = [] # Index that hold at least one frequency
t0 = [] # Initial time (ps)
frequencies = [] # Frequency (GHz)
quality = [] # Quality factor
chi = [] # Chi Squared
meanqdiff = [] # Mean Squared Difference
nterms = [] # Number of fit terms
        
# Now, begin iteration on files
for n in names:

    print("---> File {}/{}".format(names.index(n)+1, len(names)))
    
    # Load data
    t, V, details = ivs.loadNicePumpProbe(ivs.filenameToMeasureFilename(n,
                                                                        home))
    
    # Load fit parameters
    data, header, fit_params = ivs.loadTxt(ivs.filenameToFitsFilename(n, home))
    fit_params = ivu.InstancesDict(fit_params)
    del data, header
    
    # Choose data to fit
    if fit_params.use_full_mean:
        data = np.mean(V, axis=1)
    else:
        data = np.mean(V[:, fit_params.use_experiments], axis=1)

    # Make a vertical shift
    data = data - fit_params.voltage_zero

    # Choose time interval to fit
    t0n = fit_params.time_range[0] # Initial time assumed to optimize it
    i = np.argmin(np.abs(t-t0n)) # We'll take this index as main initial time

    # For each file, we'll have a different set of data to collect
    f_jmean = []
    f_jgood = [] # From here on, this is data I wouldn't like to overwrite
    f_jreallygood = []
    f_t0 = []
    f_frequencies = []
    f_quality = []
    f_chi = []
    f_meanqdiff = []
    f_nterms = []

    # Now we can iterate over the initial time
    if i-Ni < 0:
        posiblej = list(range(0, Ni))
    else:
        posiblej = list(range(i-Ni, i+Ni))
    for j in posiblej:
    
        print("Initial Time {}/{}".format(posiblej.index(j)+1, 
                                             len(posiblej)))
        
        # Choose initial time t0
        t0j = t[j]
        f_t0.append(t0j)
        
        # Crop data accorddingly
        tj, dataj = iva.cropData(t0j, t, data)
        fit_params.time_range = (t0j, t[-1])
        fit_params.voltage_zero = 0
        
        # Use linear prediction, if allowed
        try:
            results, others, plots = iva.linearPrediction(
                    tj, 
                    dataj, 
                    details['dt'], 
                    svalues=fit_params.Nsingular_values,
                    printing=False)
            f_jgood.append(j)
            fit_terms = plots.fit
            del plots

            # Keep only the fits that satisfy us
            if results.shape[0]!=1: # Select closest frequency to desired one
                imax = np.argmin(np.abs(results[:,0] - 
                                        desired_frequency * 
                                        np.ones(len(results[:,0]))))
                if results[imax,0] != 0:
                    f_frequencies.append(results[imax,0])
                    f_quality.append(results[imax,2])
                    f_chi.append(others['chi_squared'])
                    f_jreallygood.append(j)
                    f_meanqdiff.append( np.mean( (fit_terms[:,1] - 
                                                  fit_terms[:,imax+2])**2 ) )
                    f_nterms.append(results.shape[0])
            else:
                if results[0,0] != 0:
                    f_frequencies.append(results[0,0])
                    f_quality.append(results[0,2])
                    f_chi.append(others['chi_squared'])
                    f_jreallygood.append(j)
                    f_meanqdiff.append( np.mean( (fit_terms[:,1] - 
                                                  fit_terms[:,imax+2])**2 ) )
                    f_nterms.append(results.shape[0])
            
        except:
            pass
            
    del t0j, tj, dataj, posiblej
    del results, others, t, V, data, details, fit_terms

    # Now, before going to the next file, save data
    jmean.append(imax)
    jgood.append(f_jgood)
    jreallygood.append(f_jreallygood)
    t0.append(f_t0)
    frequencies.append(f_frequencies)
    quality.append(f_quality)
    chi.append(f_chi)
    meanqdiff.append(f_meanqdiff)
    nterms.append(f_nterms)

del f_jgood, f_jreallygood, f_frequencies, f_quality, f_chi, f_meanqdiff, 
del f_nterms, j, imax, n

#%%

if len(names)==1:
    jmean = jmean[0]
    jgood = jgood[0]
    jreallygood = jreallygood[0]
    frequencies = frequencies[0]
    quality = quality[0]
    chi = chi[0]
    meanqdiff = meanqdiff[0]
    nterms = nterms[0]

    # Make plots showing results
    plt.figure()
    
    plt.subplot(2,2,1)
    plt.plot(jreallygood, frequencies, 'x')
    plt.plot(jmean, frequencies[jmean], 'xr')
    plt.ylabel('Frecuencia (GHz)')
    
    ax = plt.subplot(2,2,2)
    plt.plot(jreallygood, quality, 'o')
    plt.plot(jmean, quality[jmean], 'or')
    ax.axes.yaxis.tick_right()
    plt.ylabel('Factor de calidad')
    #ax.axes.yaxis.label_position = 'right'
    
    plt.subplot(2,2,3)
    plt.plot(jreallygood, chi, '.')
    plt.plot(jmean, chi[jmean], 'xr')
    plt.ylabel('Chi cuadrado')
    
    ax = plt.subplot(2,2,4)
    plt.plot(jreallygood, meanqdiff, 'x')
    plt.plot(jmean, meanqdiff[jmean], 'xr')
    ax.axes.yaxis.tick_right()
    plt.ylabel('Diferencia cuadrática media')
    #ax.axes.yaxis.label_position = 'right'

#%%

# Save data
data = np.array([jreallygood, list(t[jreallygood]), frequencies, quality, chi, meanqdiff]).T#, stdqdiff]).T
header = ['Índice temporal inicial', 'Tiempo inicial (ps)', 'Frecuencia (GHz)', 
          'Factor de calidad', 'Chi cuadrado', 'Diferencia cuadrática media']#, 
#          'Desviación estándar de la diferencia cuadrática']
fit_params.update(dict(svalues=4, i=i, Ni=Ni))
ivs.saveTxt(os.path.join(home, r'Análisis/{}_LP.txt'.format(name)), data, 
            header=header, footer=fit_params.__dict__)