# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import iv_utilities_module as ivu
import os
import numpy as np

#%% PARAMETERS -------------------------------------------------------------------

# Parameters
name = 'M_20190408_04'
path = r'C:\Users\Luciana\Desktop\Vale e Iván\Mediciones\2019-04-08'

# Plot parameters
plot_params = dict(
        plot = False,
        interactive = False,
        autoclose = True,
        autosave = True
        )
plot_params = ivu.InstancesDict(plot_params)

# Fit parameters
fit_params = dict(
        round_Matlab_needed = True, # Pyhon 3.6.2 needs it
        use_full_mean = True,
        use_experiments = [0], # First is 0, not 1!
        send_tail_to_zero = False,
        use_fraction = .2,
        choose_t0 = True,
        choose_tf = False
        )
fit_params = ivu.InstancesDict(fit_params)

# Create full filename
filename = os.path.join(path, name+'.txt')

#%% PLOT -------------------------------------------------------------------------

# Plot
if plot_params.plot:
    ivp.plotPumpProbe(filename, 
                      interactive=plot_params.interactive, 
                      autosave=plot_params.autosave)

# Several plots
#ivp.plotAllPumpProbe(path, autosave=autosave, autoclose=autoclose)

#%% LINEAR PREDICTION -------------------------------------------------------------

# Load data
t, V, details = ivs.loadNicePumpProbe(filename)

# Choose time interval to fit
if fit_params.choose_t0: # Choose initial time t0
    t0 = ivp.interactiveTimeSelector(filename, autoclose=plot_params.autoclose)
    t, V = iva.cropData(t0, t, V)
else:
    t0 = t[0]
if fit_params.choose_tf: # Choose final time tf
    tf = ivp.interactiveTimeSelector(filename, autoclose=plot_params.autoclose)
    t, V = iva.cropData(tf, t, V, logic='<=')
else:
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
    V0 = min(data[int( (1-fit_params.use_fraction) * len(data)):]) 
else:
    V0 = 0
data = data - V0
fit_params.voltage_zero = V0
del V0

# Use linear prediction
results, other_results, plot_results = iva.linearPrediction(
    t, data, details['dt'], 
    autoclose=plot_params.autoclose,
    round_Matlab_needed=fit_params.round_Matlab_needed)

# Plot linear prediction
ivp.linearPredictionPlot(filename, plot_results, autosave=plot_params.autosave)

# Generate fit tables
tables = iva.linearPredictionTables(fit_params, results, other_results)
ivu.copy(tables[0])