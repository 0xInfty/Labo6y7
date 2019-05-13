# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import os
import numpy as np

#%% PARAMETERS -------------------------------------------------------------------

# Parameters
name = 'M_20190508_02'
path = r'C:\Users\Luciana\Desktop\Vale e Iván\Mediciones\2019-05-08'

# Plot parameters
plot = True
interactive = False
autoclose = True
autosave = True

# Fit parameters
round_Matlab_needed = True # Pyhon 3.6.2 needs it
use_full_mean = False
use_experiments = [2] # First is 0, not 1!
send_tail_to_zero = True
use_fraction = .2
choose_tf = True

# Create full filename
filename = os.path.join(path, name+'.txt')

#%% PLOT -------------------------------------------------------------------------

# Plot
if plot:
    ivp.plotPumpProbe(filename, interactive=interactive, autosave=autosave)

# Several plots
#ivp.plotAllPumpProbe(path, autosave=autosave, autoclose=autoclose)

#%% LINEAR PREDICTION -------------------------------------------------------------

# Load data
t, V, details = ivs.loadNicePumpProbe(filename)
t0 = ivp.interactiveTimeZero(filename, autoclose=autoclose)
t, V = iva.cropData(t0, t, V)
if choose_tf:
    tf = ivp.interactiveTimeZero(filename, autoclose)
    t, V = iva.cropData(tf, t, V, logic='<=')
dt = details['dt']

# Use linear prediction
if use_full_mean:
    data = np.mean(V, axis=1)
else:
    data = np.mean(V[:, use_experiments], axis=1)
V0 = min(data[int( (1-use_fraction) * len(data)):])
data = data - V0
results, others = iva.linearPrediction(t, data, dt, 
                                      autoclose=autoclose,
                                      round_Matlab_needed=round_Matlab_needed)

# Plot linear prediction
fig = ivp.linearPredictionPlot(filename, others, autosave=autosave)
