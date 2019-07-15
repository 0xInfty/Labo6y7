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
home = r'C:\Users\Vall\OneDrive\Labo 6 y 7'
names = ['M_20190605_11']
desired_frequency = 9 # Desired frequency for the ideal fit
Ni = 40 # How many index around the main one we'll try for the initial time

#%%

# Data to collect while iterating
jmean = [] # Mean index
jgood = [] # Index that allow fitting
jreallygood = [] # Index that hold at least one frequency
t0 = [] # Initial time (ps)
data0 = []
t = []
data = []
frequencies = [] # Frequency (GHz)
quality = [] # Quality factor
chi = [] # Chi Squared
meanqdiff = [] # Mean Squared Difference
nterms = [] # Number of fit terms
        
# Now, begin iteration on files
for n in names:

    print("---> File {}/{}".format(names.index(n)+1, len(names)))
    
    # Load data
    t_n, V, details = ivs.loadNicePumpProbe(
        ivs.filenameToMeasureFilename(n,home))
    
    # Load fit parameters
    results, header, fit_params = ivs.loadTxt(
        ivs.filenameToFitsFilename(n, home))
    fit_params = ivu.InstancesDict(fit_params)
    del results, header
    
    # Choose data to fit
    if fit_params.use_full_mean:
        data_n = np.mean(V, axis=1)
    else:
        data_n = np.mean(V[:, fit_params.use_experiments], axis=1)

    # Make a vertical shift
    data_n = data_n - fit_params.voltage_zero

    # Choose time interval to fit
    t0_n = fit_params.time_range[0] # Initial time assumed to optimize it
    i = np.argmin(np.abs(t_n-t0_n)) # We'll take this index as main initial time

    # For each file, we'll have a different set of data to collect
    jgood_n = [] # From here on, this is data I wouldn't like to overwrite
    jreallygood_n = []
    t0_n = []
    frequencies_n = []
    quality_n = []
    chi_n = []
    meanqdiff_n = []
    nterms_n = []

    # Now we can iterate over the initial time
    if i-Ni < 0:
        posiblej = list(range(0, Ni))
    else:
        posiblej = list(range(i-Ni, i+Ni))
    t0.append(t_n[posiblej])
    data0.append(data_n[posiblej])
    for j in posiblej:
    
        print("Initial Time {}/{}".format(posiblej.index(j)+1, 
                                             len(posiblej)))
        
        # Choose initial time t0
        t0_j = t_n[j]
        t0_n.append(t0_j)
        
        # Crop data accorddingly
        t_j, data_j = iva.cropData(t0_j, t_n, data_n)
        fit_params.time_range = (t_j[0], t_j[-1])
        fit_params.voltage_zero = 0
        
        # Use linear prediction, if allowed
        try:
            results, others, plots = iva.linearPrediction(
                    t_j, 
                    data_j,
                    details['dt'], 
                    svalues=fit_params.Nsingular_values,
                    printing=False)
            jgood_n.append(j)
            fit_terms = plots.fit
            del plots

            # Keep only the fits that satisfy us
            if results.shape[0]!=1: # Select closest frequency to desired one
                imax = np.argmin(np.abs(results[:,0] - 
                                        desired_frequency * 
                                        np.ones(len(results[:,0]))))
                if results[imax,0] != 0:
                    frequencies_n.append(results[imax,0])
                    quality_n.append(results[imax,2])
                    chi_n.append(others['chi_squared'])
                    jreallygood_n.append(j)
                    meanqdiff_n.append( np.mean( (fit_terms[:,1] - 
                                                  fit_terms[:,imax+2])**2 ) )
                    nterms_n.append(results.shape[0])
            else:
                if results[0,0] != 0:
                    frequencies_n.append(results[0,0])
                    quality_n.append(results[0,2])
                    chi_n.append(others['chi_squared'])
                    jreallygood_n.append(j)
                    meanqdiff_n.append( np.mean( (fit_terms[:,1] - 
                                                  fit_terms[:,imax+2])**2 ) )
                    nterms_n.append(results.shape[0])
            
        except:
            pass
            
    del j, t0_j, t_j, data_j, posiblej
    del results, others, V, details, fit_terms

    # Now, before going to the next file, save data
    jmean.append(i)
    jgood.append(jgood_n)
    jreallygood.append(jreallygood_n)
    t0.append(t0_n)
    t.append(t_n)
    data.append(data_n)
    frequencies.append(frequencies_n)
    quality.append(quality_n)
    chi.append(chi_n)
    meanqdiff.append(meanqdiff_n)
    nterms.append(nterms_n)

del jgood_n, jreallygood_n, t_n, data_n, t0_n
del frequencies_n, quality_n, chi_n, meanqdiff_n, nterms_n
del i, imax, n

#%%

if len(names)==1:
    jmean = jmean[0]
    jgood = jgood[0]
    jreallygood = jreallygood[0]
    t = t[0]
    data = data[0]
    t0 = t0[0]
    data0 = data0[0]
    frequencies = frequencies[0]
    quality = quality[0]
    chi = chi[0]
    meanqdiff = meanqdiff[0]
    nterms = nterms[0]
    names = names[0]

#%%

    # Make a general plot showing the chosen initial times
    plt.figure()
    ax = plt.subplot()
    plt.plot(t, data, 'k', linewidth=0.5)
    plt.plot(t0, data0, 'r')
    plt.ylabel(r'Voltaje ($\mu$V)')
    plt.xlabel(r'Tiempo (ps)')
    ax.minorticks_on()
    ax.tick_params(axis='y', which='minor', left=False)
    ax.tick_params(length=5)
    ax.grid(axis='x', which='both')

    # Make plots showing results
    fig = plt.figure()
    grid = plt.GridSpec(5, 1, hspace=0)
    
    # Voltage plot
    ax0 = plt.subplot(grid[0,0])
    plt.plot(t0, data0, 'k')
    ax0.axes.xaxis.tick_top()
    ax0.minorticks_on()
    ax0.tick_params(axis='y', which='minor', length=0)
    ax0.tick_params(length=5)
    ax0.set_xlabel('Tiempo inicial (ps)')
    ax0.axes.xaxis.set_label_position('top')
    ax0.set_ylabel(r'Voltaje ($\mu$s)')
    ax0.grid(axis='x', which='both')
    plt.show()
    xlim = ax0.get_xlim()
    
    # Frequency plot, right axis
    ax1 = plt.subplot(grid[1:4,0])
    plt.plot(np.array(t0)[jreallygood], frequencies, 'or')
    ax1.set_xlim(xlim)
    ax1.axes.xaxis.tick_top()
    ax1.minorticks_on()
    ax1.set_ylabel('Frecuencia (GHz)', color='tab:red')
    ax1.tick_params(axis='y', labelcolor='tab:red')
    ax1.tick_params(axis='y', which='minor', length=0)
    ax1.grid(axis='x', which='both')
    
    # Quality factor, left axis
    ax2 = ax1.twinx()  # Second axes that shares the same x-axis
    ax2.set_ylabel('Factor de calidad (u.a.)', color='tab:blue')
    plt.plot(np.array(t0)[jreallygood], quality, 'xb', markersize=7)
    ax2.tick_params(axis='y', labelcolor='tab:blue')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    for l in ax1.get_xticklabels():
        l.set_visible(False)
    del l
    
    # Number of terms
    ax3 = plt.subplot(grid[-1,0])
    plt.plot(np.array(t0)[jreallygood], nterms, 'og')
    ax3.set_xlim(xlim)
    ax3.minorticks_on()
    ax3.tick_params(axis='y', which='minor', left=False)
    ax3.tick_params(length=5)
    ax3.grid(axis='x', which='both')
    for l in ax3.get_xticklabels():
        l.set_visible(False)
    del l
    ax3.set_ylabel("Número de \ntérminos")
    
    # Mean initial time
    ylim = ax0.get_ylim()
    ax0.vlines(t0[jmean], ylim[0], ylim[1], linewidth=1)
    ax0.set_ylim(ylim)
    ylim = ax1.get_ylim()
    ax1.vlines(t0[jmean], ylim[0], ylim[1], linewidth=1)
    ax1.set_ylim(ylim)
    ylim = ax3.get_ylim()
    ax3.vlines(t0[jmean], ylim[0], ylim[1], linewidth=1)
    ax3.set_ylim(ylim)
    del ylim
    
    # Make plots showing statistics
    fig = plt.figure()
    grid = plt.GridSpec(5, 1, hspace=0)
    
    # Voltage plot
    ax0 = plt.subplot(grid[0,0])
    plt.plot(t0, data0, 'k')
    ax0.axes.xaxis.tick_top()
    ax0.minorticks_on()
    ax0.tick_params(axis='y', which='minor', length=0)
    ax0.tick_params(length=5)
    ax0.set_xlabel('Tiempo inicial (ps)')
    ax0.axes.xaxis.set_label_position('top')
    ax0.set_ylabel(r'Voltaje ($\mu$s)')
    ax0.grid(axis='x', which='both')
    plt.show()
    xlim = ax0.get_xlim()
    
    # Chi Squared
    ax1 = plt.subplot(grid[1:3,0])
    plt.plot(np.array(t0)[jreallygood], chi, 'or')
    ax1.set_xlim(xlim)
#    ax1.axes.yaxis.label_position = 'right'
    ax1.axes.yaxis.tick_right()
    ax1.minorticks_on()
    ax1.set_ylabel('Chi cuadrado')
    ax1.tick_params(axis='y')
    ax1.tick_params(axis='y', which='minor', length=0)
    ax1.grid(axis='x', which='both')
    
    # Mean Squared Difference
    ax2 = plt.subplot(grid[3:,0])
    plt.plot(np.array(t0)[jreallygood], meanqdiff, 'ob')
    ax2.set_xlim(xlim)
    ax2.minorticks_on()
    ax2.set_ylabel('Diferencia \ncuadrática media')
    ax2.tick_params(axis='y')
    ax2.tick_params(axis='y', which='minor', length=0)
    ax2.grid(axis='x', which='both')
    plt.show()
    for l in ax1.get_xticklabels():
        l.set_visible(False)
    del l
    
    # Mean initial time
    ylim = ax0.get_ylim()
    ax0.vlines(t0[jmean], ylim[0], ylim[1], linewidth=1)
    ax0.set_ylim(ylim)
    ylim = ax1.get_ylim()
    ax1.vlines(t0[jmean], ylim[0], ylim[1], linewidth=1)
    ax1.set_ylim(ylim)
    ylim = ax2.get_ylim()
    ax2.vlines(t0[jmean], ylim[0], ylim[1], linewidth=1)
    ax2.set_ylim(ylim)
    del ylim

#%%

    # Save data
    data = np.array([jreallygood, list(t[jreallygood]), 
                     frequencies, quality, chi, meanqdiff]).T#, stdqdiff]).T
    header = ['Índice temporal inicial', 'Tiempo inicial (ps)', 'Frecuencia (GHz)', 
              'Factor de calidad', 'Chi cuadrado', 'Diferencia cuadrática media']#, 
    #          'Desviación estándar de la diferencia cuadrática']
    fit_params.update(dict(svalues=4, i=jmean, Ni=Ni))
    ivs.saveTxt(os.path.join(home, r'Análisis/{}_LP.txt'.format(names)), data, 
                header=header, footer=fit_params.__dict__)