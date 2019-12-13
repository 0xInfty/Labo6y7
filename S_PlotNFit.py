# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import iv_utilities_module as ivu
import numpy as np

#%% PARAMETERS -------------------------------------------------------------------

# Parameters
name = 'M_20191129_01'
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'

# Save parameters
autosave = True
overwrite = True

# Plot parameters
plot_params = dict(
        plot = False,
        interactive = True,
        autoclose = True,
        extension = '.png'
        )
plot_params = ivu.InstancesDict(plot_params)

# Fit parameters
fit_params = dict(
        use_full_mean = True,
        use_experiments = [1], # First is 0, not 1!
        send_tail_to_zero = True,
        tail_method = 'mean', # Could also be 'min' or 'max' or any numpy function
        use_fraction = .1,
        choose_t0 = True,
        choose_tf = False,
        svalues = None,
        max_svalues = 20,
        )
fit_params = ivu.InstancesDict(fit_params)

# Create full filename
filename = ivs.filenameToMeasureFilename(name, home=home)

#%% PLOT --------------------------------------------------------------------------

# Plot
if plot_params.plot:
    fig, legb, savb = ivp.plotPumpProbe(filename,
                                        interactive=plot_params.interactive, 
                                        extension=plot_params.extension,
                                        autosave=False,
                                        overwrite=True,
#                                        loc='upper right'
                                        )

if False: print(
""" TO PLOT SEVERAL MEASUREMENTS
import os
path = os.path.split(filename)[0]
ivp.plotAllPumpProbe(path,
                     autoclose=plot_params.autoclose,
                     extension=plot_params.extension,
                     autosave=autosave,
                     overwrite=True)
"""
    )

#%% LINEAR PREDICTION -------------------------------------------------------------

# Load data
t, V, details = ivs.loadNicePumpProbe(filename)

# Choose time interval to fit
if fit_params.choose_t0: # Choose initial time t0
    t0 = ivp.interactiveTimeSelector(filename, autoclose=plot_params.autoclose)
    t, V = iva.cropData(t0, t, V)
else:
    try:
        t, V = iva.cropData(t0, t, V)
    except NameError:
        t0 = t[0]
if fit_params.choose_tf: # Choose final time tf
    tf = ivp.interactiveTimeSelector(filename, autoclose=plot_params.autoclose)
    t, V = iva.cropData(tf, t, V, logic='<=')
else:
    try:
        t, V = iva.cropData(tf, t, V, logic='<=')
    except NameError:
        tf = t[-1]
fit_params.time_range = (t0, tf)
del t0, tf

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
del V0

# Use linear prediction
results, other_results, plot_results = iva.linearPrediction(
    t, data, details['dt'],
    svalues=fit_params.svalues,
    max_svalues=fit_params.max_svalues,
    autoclose=plot_params.autoclose)
if autosave:
    ivs.linearPredictionSave(filename, results, other_results, fit_params,
                             overwrite=overwrite)
 
# Plot linear prediction
ivp.linearPredictionPlot(filename, plot_results, 
                         autosave=autosave,
                         extension=plot_params.extension,
                         overwrite=overwrite)

# Generate fit tables
tables = iva.linearPredictionTables(fit_params, results, other_results)
ivu.copy(tables[0])