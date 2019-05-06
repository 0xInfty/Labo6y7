# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import os
#import numpy as np

#%% CHOOSE FILE -------------------------------------------------------------------

# Parameters
name = 'M_20190506_11'
folder = '2019-05-06'
path = r'F:\Pump-Probe\Iván y Valeria\Mediciones'

# Create full filename
filename = os.path.join(path, folder, name+'.txt')

#%% ONE PLOT ----------------------------------------------------------------------

# Parameters
full = False
autosave = True

# Plot
if full:
    ivp.plotFullPumpProbe(filename, autosave=autosave)
else:
    ivp.plotPumpProbe(filename, autosave=autosave)

#%% SEVERAL PLOTS -----------------------------------------------------------------

# Parameters
full = True
autosave = True
autoclose = False

# Plot
ivp.plotAllPumpProbe(path, full=full, 
                     autosave=autosave, autoclose=autoclose)

#%% LINEAR PREDICTION -------------------------------------------------------------

# Parameters
autoclose = True
autosave = True
round_Matlab_needed = True # Pyhon 3.6.2 needs it

# Load data
t, V, meanV, details = ivs.loadNicePumpProbe(filename)
t0 = ivp.interactiveTimeZero(filename, autoclose=autoclose)
t, V, meanV = iva.cropData(t0, t, V, meanV)
#t, V, meanV, details = iva.loadZeroPumpProbe(filename)
dt = details['dt']

# Use linear prediction
data = meanV
results, others = iva.linearPrediction(t, data, dt, 
                                      autoclose=autoclose,
                                      round_Matlab_needed=round_Matlab_needed)

# Plot linear prediction
fig = ivp.linearPredictionPlot(filename, others, autosave=autosave)
