# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 09:48:07 2019

@author: Valeria
"""

import iv_save_module as ivs
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
import numpy as np
import os

#%% LINEAR PREDICTION REPRISE

# Parameters
names = ['M_20190612_02', 'M_20191025_05']
labels = ['Nano-paralelepípedo', 'Nano-cilindro']
scale = [2, 1]
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'

# Save parameters
path = os.path.join(home)
overwrite = False
autosave = True

# Plot parameters
plot_params = dict(
        plot = True,
        interactive = False,
        autoclose = True,
        extension = '.png'
        )
plot_params = ivu.InstancesDict(plot_params)

#%% MORE PARAMETERS

# Define other parameters
filenames = [ivs.filenameToMeasureFilename(n, home=home) for n in names]
fit_filenames = [ivs.filenameToFitsFilename(n, home=home) for n in names]

#%% LOAD DATA

# Load data from files
data = []
for f in filenames:
    t, V, details = ivs.loadNicePumpProbe(f)
    data.append(np.array([t,*V.T]).T)
del t, V, details, f

# Load data from fit filenames
other_results_keys = ['Nsingular_values', 'chi_squared']
fit_params = []
for f in fit_filenames:
    
    # Load data from a base fit made by hand
    r, fit_header, ft = ivs.loadTxt(f)
    
    # Reorganize data
    others = {k: ft[k] for k in other_results_keys}
    fp = dict(ft)
    for k in other_results_keys:
        fp.pop(k)
    del k
    fp = ivu.InstancesDict(fp)
    del ft
    
    # Add data to external variables
    fit_params.append(fp)

del r, others, fp, f

#%% MAKE PLOT

# Use linear prediction

results = []
other_results = []

r, ores, pr = iva.linearPrediction(
    t, data, details['dt'],
    svalues=fit_params.svalues,
    max_svalues=fit_params.max_svalues,
    autoclose=plot_params.autoclose)

#%%

#if len(names) > 2:
#    raise ValueError("No está armado para más de 2 mediciones :P")
#
#plt.figure()
#grid = plt.GridSpec(2, 1, hspace=0.2)
#axp, axf = [plt.subplot(g) for g in grid]
#
#axp = [axp, axp.twinx()]
#colors = ['r', 'b']
#
#for a, d, fd, fp, c, l in zip(axp, data, results, fit_params, colors, labels):
#    
#    # Plot data
#    t = d[:,0]
#    Vmean = np.mean(d[:,1:], axis=1)
#    a.plot(t, Vmean, color=c, linewidth=.3)
#    a.tick_params(axis='y', labelcolor=c, )
#    a.set_ylabel(r'Voltaje ($\mu$V)', color=c)
#    
#    # Plot fit
#    t0 =  fp.time_range[0]
#    try:
#        omega = 2 * np.pi * fd[:,0]
#        tau = fd[:,1]
#        damping_time = 1/tau
#        amplitudes = fd[:,3]
#        phases = fd[:,4]
#    except:
#        omega = np.array([2 * np.pi * fd[0]])
#        tau = np.array([fd[1]])
#        damping_time = 1/tau
#        amplitudes = np.array([fd[3]])
#        phases = np.array([fd[4]])
#    fit_terms = np.array([amp * np.exp(-b*(t-t0)) * np.cos(w*(t-t0) + phi)
#                         for amp, b, w, phi in zip(amplitudes,
#                                                   damping_time,
#                                                   omega,
#                                                   phases)]).T
#    fit = sum(fit_terms.T)
#    a.plot(t, fit, color=c, linewidth=1.5)