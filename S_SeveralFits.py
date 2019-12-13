# -*- coding: utf-8 -*-
"""
Created on Mon Dec 09 14:28:53 2019

@author: Vall
"""

from itertools import combinations as comb
import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import iv_utilities_module as ivu
from matplotlib.pyplot import close
import numpy as np
import os

#%% LINEAR PREDICTION REPRISE

# Parameters
name = 'M_20191018_11'
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
nexperiments = 4
groups_mode = 'own' # Combinations 'comb', each experiment by its own 'own'
series = 'Power'

# Save parameters
path = os.path.join(home, 'Análisis', series + '_' + name)
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

#%% LOAD DATA

# Make filenames routs
filename = ivs.filenameToMeasureFilename(name, home=home)
fit_filename = ivs.filenameToFitsFilename(name, home=home)
other_fit_filename = os.path.join(path, name+'.txt')

# Load data from a base fit made by hand
results, header, footer = ivs.loadTxt(fit_filename)

# Reorganize data
other_results_keys = ['Nsingular_values', 'chi_squared']
other_results = {k: footer[k] for k in other_results_keys}
fit_params = dict(footer)
for k in other_results_keys:
    fit_params.pop(k)
fit_params = ivu.InstancesDict(fit_params)
del footer

# New parameters
fit_params.use_full_mean = False
fit_params.choose_t0 = False
fit_params.choose_tf = False
t0 = fit_params.time_range[0]
tf = fit_params.time_range[-1]
fit_params.svalues = other_results['Nsingular_values']

#%% MAKE SEVERAL FITS

# Make groups of experiments
if groups_mode == 'comb':
    # Make all combinations of 3 elements
    experiments_groups = [list(c) for c in comb(list(range(nexperiments)), 3)]
elif groups_mode == 'own':
    # Take each experiment on its own
    experiments_groups = [[i] for i in range(nexperiments)]
print("Watch it! Loop will make " + str(len(experiments_groups)) +
      " repetitions. Press Ctrl+C to cancel.")

all_failed_groups = []
all_groups = []
all_results = []
all_other_results = []
all_fit_params = []
#all_tables = []
for i, g in enumerate(experiments_groups):

    print("----> Loop {}/{}. Experiments: {}".format(
            i+1, len(experiments_groups), g))
    
    try:
    
        fit_params.use_experiments = g
        
        # Load data
        t, V, details = ivs.loadNicePumpProbe(filename)
        
        # Choose time interval to fit
        if fit_params.choose_t0: # Choose initial time t0
            t0 = ivp.interactiveTimeSelector(filename, 
                                             autoclose=plot_params.autoclose)
            t, V = iva.cropData(t0, t, V)
        else:
            try:
                t, V = iva.cropData(t0, t, V)
            except NameError:
                t0 = t[0]
        if fit_params.choose_tf: # Choose final time tf
            tf = ivp.interactiveTimeSelector(filename, 
                                             autoclose=plot_params.autoclose)
            t, V = iva.cropData(tf, t, V, logic='<=')
        else:
            try:
                t, V = iva.cropData(tf, t, V, logic='<=')
            except NameError:
                tf = t[-1]
        fit_params.time_range = (t0, tf)
        
        # Choose data to fit
        if fit_params.use_full_mean:
            data = np.mean(V, axis=1)
        else:
            data = np.mean(V[:, fit_params.use_experiments], axis=1)
        
        # Make a vertical shift
        if fit_params.send_tail_to_zero:
            function = eval('np.{}'.format(fit_params.tail_method))
            V0 = function(data[int( (1-fit_params.use_fraction) * len(data)):])
            del function
        else:
            try:
                V0
            except NameError:
                V0 = 0
        data = data - V0
        fit_params.voltage_zero = V0
        
        # Use linear prediction
        results, other_results, plot_results = iva.linearPrediction(
            t, data, details['dt'],
            svalues=fit_params.svalues,
            autoclose=plot_params.autoclose)
        if autosave:
            ivs.linearPredictionSave(other_fit_filename, results, other_results, 
                                     fit_params, overwrite=overwrite)
         
        # Plot linear prediction
        ivp.linearPredictionPlot(other_fit_filename, plot_results, 
                                 autosave=autosave,
                                 extension=plot_params.extension,
                                 overwrite=overwrite)
        close()
        
#        # Generate fit tables
#        tables = iva.linearPredictionTables(fit_params, 
#                                            results, 
#                                            other_results))
        
        # Save data
        all_groups.append(g)
        all_results.append(results)
        all_other_results.append(other_results)
        all_fit_params.append(fit_params)
#        all_tables.append(tables)
        
    except:
        
        all_failed_groups.append(g)
        print("Failed! :(")

del i, g, experiments_groups, t0, tf, V0, t, data
del results, other_results, plot_results, fit_params, plot_params