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
import random as ran

#%%

# Parameters
home = r'C:\Users\Vall\OneDrive\Labo 6 y 7'
rods_filename = os.path.join(home, r'Análisis\Rods_LIGO1.txt')
#sem_filename = os.path.join(home, r'Muestras\SEM\LIGO1\LIGO1 Geometrías\1\Resultados_SEM_LIGO1_1.txt')
desired_frequency = 9 # Desired frequency for the ideal fit
Ni = 40 # How many index around the main one we'll try for the initial time
autosave = True
autoclose = True

## --> Rare Series
#names = ['M_20190610_07', 'M_20190605_07', 'M_20190610_13', 'M_20190610_01', 'M_20190610_12'] # OUTLIERS
#series = 'Rare'

# --> Random Series
#names = ['M_20190605_07', 'M_20190605_11', 'M_20190605_12', 'M_20190610_06', 'M_20190610_07']
#series = 'Random_1'

# Look for the list of rods and filenames
filenames = [] # Will contain filenames like 'M_20190610_01'
rods = [] # Will contain rods' positions like '1,2'
with open(rods_filename, 'r') as file:
    for line in file:
        if line[0]!='#':
            filenames.append(line.split('\t')[0]) # Save filenames
            rods.append(line.split('\t')[1].split('\n')[0]) # Save rods
    del line

index = ran.sample(range(len(rods)), 5)
names = [filenames[i] for i in index]
series = 'Random_1'

# Keep only the selected filenames and rods
index = [filenames.index(n) for n in names]
rods = [rods[i] for i in index]
del filenames

#%%

def filenameToFigFilename(filename, series='', home=home):
    
    """Given a filename 'M_20190610_01', returns path to fits' data"""    

    if series!='':
        series = '_{}'.format(series)
    base = os.path.join(home, r'Análisis/StudyLP'+series)
    if not os.path.isdir(base):
        os.makedirs(base)

    date = filename.split('_')[1] # From 'M_20190610_01' take '20190610'
    date = '-'.join([date[:4], date[4:6], date[6:]]) # Transfrom to '2019-06-10'
            
    fig_filenames = [
            os.path.join(base, filename+'_Voltage.png'),
            os.path.join(base, filename+'_Params.png'),
            os.path.join(base, filename+'_Stats.png')
            ]
    
    return fig_filenames

def filenameToFilename(filename, series='', home=home):
    
    """Given a filename 'M_20190610_01', returns path to fits' data"""
    
    if series!='':
        series = '_{}'.format(series)
    base = os.path.join(home, r'Análisis/StudyLP'+series)
    if not os.path.isdir(base):
        os.makedirs(base)
    
    date = filename.split('_')[1] # From 'M_20190610_01' take '20190610'
    date = '-'.join([date[:4], date[4:6], date[6:]]) # Transfrom to '2019-06-10'
    filename = os.path.join(base, filename+'.txt')
    
    return filename

def figsFilename(fig_name, series=''):
    
    """Given a fig_name 'DifCuadráticaMedia', returns path to fig"""
    
    if series!='':
        series = '_{}'.format(series)
    base = os.path.join(home, r'Análisis/StudyLP'+series)
    if not os.path.isdir(base):
        os.makedirs(base)
    
    filename = os.path.join(base, fig_name+'.png')
    
    return filename

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
fit_params = []
        
# Now, begin iteration on files
for n in names:

    print("---> File {}/{}".format(names.index(n)+1, len(names)))
    
    # Load data
    t_n, V, details = ivs.loadNicePumpProbe(
        ivs.filenameToMeasureFilename(n,home))
    
    # Load fit parameters
    results, header, fit_params_n = ivs.loadTxt(
        ivs.filenameToFitsFilename(n, home))
    fit_params_n = ivu.InstancesDict(fit_params_n)
    del results, header
    
    # Choose data to fit
    if fit_params_n.use_full_mean:
        data_n = np.mean(V, axis=1)
    else:
        data_n = np.mean(V[:, fit_params_n.use_experiments], axis=1)

    # Make a vertical shift
    data_n = data_n - fit_params_n.voltage_zero

    # Choose time interval to fit
    t0_n = fit_params_n.time_range[0] # Initial time assumed to optimize it
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
    if i-Ni//2 < 0:
        posiblej = list(range(0, Ni))
    else:
        posiblej = list(range(i-Ni//2, i+Ni//2))
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
        
        # Use linear prediction, if allowed
        try:
            results, others, plots = iva.linearPrediction(
                    t_j, 
                    data_j,
                    details['dt'], 
                    svalues=fit_params_n.Nsingular_values,
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
    t.append(t_n)
    data.append(data_n)
    frequencies.append(frequencies_n)
    quality.append(quality_n)
    chi.append(chi_n)
    meanqdiff.append(meanqdiff_n)
    nterms.append(nterms_n)
    fit_params.append(fit_params_n)

del jgood_n, jreallygood_n, t_n, data_n, t0_n
del frequencies_n, quality_n, chi_n, meanqdiff_n, nterms_n
del i, imax, n

#%%

for k in range(len(names)):

    # Make a general plot showing the chosen initial times
    plt.figure()
    ax = plt.subplot()
    plt.plot(t[k], data[k], 'k', linewidth=0.5)
    plt.plot(t0[k], data0[k], 'r')
    plt.ylabel(r'Voltaje ($\mu$V)')
    plt.xlabel(r'Tiempo (ps)')
    ax.minorticks_on()
    ax.tick_params(axis='y', which='minor', left=False)
    ax.tick_params(length=5)
    ax.grid(axis='x', which='both')

    # Save pitcure
    if autosave:
        plt.savefig(filenameToFigFilename(names[k], series)[0], 
                    bbox_inches='tight')
    if autoclose:
        plt.close(plt.gcf())

    # Make plots showing results
    fig = plt.figure()
    grid = plt.GridSpec(5, 1, hspace=0)
    
    # Voltage plot
    ax0 = plt.subplot(grid[0,0])
    plt.plot(t0[k], data0[k], 'k')
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
    plt.plot(t[k][jreallygood[k]], frequencies[k], 'or')
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
    plt.plot(t[k][jreallygood[k]], quality[k], 'xb', markersize=7)
    ax2.tick_params(axis='y', labelcolor='tab:blue')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    for l in ax1.get_xticklabels():
        l.set_visible(False)
    del l
    
    # Number of terms
    ax3 = plt.subplot(grid[-1,0])
    plt.plot(t[k][jreallygood[k]], nterms[k], 'og')
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
    ax0.vlines(t[k][jmean[k]], ylim[0], ylim[1], linewidth=1)
    ax0.set_ylim(ylim)
    ylim = ax1.get_ylim()
    ax1.vlines(t[k][jmean[k]], ylim[0], ylim[1], linewidth=1)
    ax1.set_ylim(ylim)
    ylim = ax3.get_ylim()
    ax3.vlines(t[k][jmean[k]], ylim[0], ylim[1], linewidth=1)
    ax3.set_ylim(ylim)
    del ylim
    
    # Save pitcure
    if autosave:
        plt.savefig(filenameToFigFilename(names[k], series)[1], 
                    bbox_inches='tight')
    if autoclose:
        plt.close(plt.gcf())

    # Make plots showing statistics
    fig = plt.figure()
    grid = plt.GridSpec(5, 1, hspace=0)
    
    # Voltage plot
    ax0 = plt.subplot(grid[0,0])
    plt.plot(t0[k], data0[k], 'k')
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
    plt.plot(t[k][jreallygood[k]], chi[k], 'or')
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
    plt.plot(t[k][jreallygood[k]], meanqdiff[k], 'ob')
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
    ax0.vlines(t[k][jmean[k]], ylim[0], ylim[1], linewidth=1)
    ax0.set_ylim(ylim)
    ylim = ax1.get_ylim()
    ax1.vlines(t[k][jmean[k]], ylim[0], ylim[1], linewidth=1)
    ax1.set_ylim(ylim)
    ylim = ax2.get_ylim()
    ax2.vlines(t[k][jmean[k]], ylim[0], ylim[1], linewidth=1)
    ax2.set_ylim(ylim)
    del ylim

    # Save pitcure
    if autosave:
        plt.savefig(filenameToFigFilename(names[k], series)[2], 
                    bbox_inches='tight')
    if autoclose:
        plt.close(plt.gcf())
        
    # Save data
    results = np.array([jreallygood[k], list(t[k][jreallygood[k]]), 
                     frequencies[k], quality[k], chi[k], meanqdiff[k]]).T#, stdqdiff]).T
    header = ['Índice temporal inicial', 'Tiempo inicial (ps)', 'Frecuencia (GHz)', 
              'Factor de calidad', 'Chi cuadrado', 'Diferencia cuadrática media']#, 
    #          'Desviación estándar de la diferencia cuadrática']
    fit_params[k].update(dict(i=jmean[k], Ni=Ni))
    ivs.saveTxt(filenameToFilename(names[k], series), results, 
                header=header, footer=fit_params[k].__dict__)

del header, results

#%% Analyse this data

# Load data
data = []
footer = []
for n in names:
    d, header, f = ivs.loadTxt(filenameToFilename(n, series))
    data.append(d)
    footer.append(f)
del d, f

# Look for the list of rods and filenames
filenames = [] # Will contain filenames like 'M_20190610_01'
rods = [] # Will contain rods' positions like '1,2'
with open(rods_filename, 'r') as file:
    for line in file:
        if line[0]!='#':
            filenames.append(line.split('\t')[0]) # Save filenames
            rods.append(line.split('\t')[1].split('\n')[0]) # Save rods
    del line

## Also load data from SEM dimension analysis
#sem_data, sem_header, sem_footer = ivs.loadTxt(sem_filename)
#other_rods = sem_footer['rods']
#new_data = []
#for r in rods:
#    i = other_rods.index(r)
#    new_data.append(sem_data[i])
#sem_data = np.array(new_data)
#del new_data, sem_footer

# Keep only data related to my selected files
index = [filenames.index(n) for n in names]
rods = [rods[i] for i in index]
#index = [other_rods.index(r) for r in rods]
#sem_data = sem_data[index,:]
del index, n, filenames

# Make several plots
plt.figure()
ax = plt.subplot()
for d in data:
    ax.plot(d[:,1]-d[0,1], d[:,2])
plt.legend(rods)
plt.xlabel('Tiempo inicial relativo (ps)')
plt.ylabel('Frecuencia (GHz)')
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='x', which='both')
plt.show()
for l in ax1.get_xticklabels():
    l.set_visible(False)
del l
if autosave:
    plt.savefig(figsFilename('Frecuencia', series), 
                bbox_inches='tight')
if autoclose:
    plt.close(plt.gcf())

plt.figure()
ax = plt.subplot()
for d in data:
    ax.plot(d[:,1]-d[0,1], d[:,3])
plt.legend(rods)
plt.xlabel('Tiempo inicial relativo (ps)')
plt.ylabel('Factor de calidad (GHz)')
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='x', which='both')
plt.show()
for l in ax1.get_xticklabels():
    l.set_visible(False)
del l
if autosave:
    plt.savefig(figsFilename('FCalidad', series), 
                bbox_inches='tight')
if autoclose:
    plt.close(plt.gcf())

plt.figure()
ax = plt.subplot()
for d in data:
    ax.plot(d[:,1]-d[0,1], d[:,4])
plt.legend(rods)
plt.xlabel('Tiempo inicial relativo (ps)')
plt.ylabel('Chi cuadrado')
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='x', which='both')
plt.show()
for l in ax1.get_xticklabels():
    l.set_visible(False)
del l
if autosave:
    plt.savefig(figsFilename('ChiCuadrado', series), 
                bbox_inches='tight')
if autoclose:
    plt.close(plt.gcf())

plt.figure()
ax = plt.subplot()
for d in data:
    ax.plot(d[:,1]-d[0,1], d[:,5])
plt.legend(rods)
plt.xlabel('Tiempo inicial relativo (ps)')
plt.ylabel('Diferencia cuadrática media')
ax.minorticks_on()
ax.tick_params(axis='y')
ax.tick_params(axis='y', which='minor', length=0)
ax.grid(axis='x', which='both')
plt.show()
for l in ax1.get_xticklabels():
    l.set_visible(False)
del l
if autosave:
    plt.savefig(figsFilename('DifCuadrática', series), 
                bbox_inches='tight')
if autoclose:
    plt.close(plt.gcf())