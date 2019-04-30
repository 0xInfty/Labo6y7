# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import os

#%% ONE PLOT ----------------------------------------------------------------------

# Parameters
name = 'M_20190408_03'
path = os.getcwd()
full = False
autosave = True

# Create full filename
filename = os.path.join(path, name+'.txt')

# Plot
if full:
    ivp.plotFullPumpProbe(filename, autosave=autosave)
else:
    ivp.plotPumpProbe(filename, autosave=autosave)

#%% SEVERAL PLOTS -----------------------------------------------------------------

# Parameters
folder = '2019-04-10'
path = os.getcwd()
full = True
autosave = True
autoclose = False

# Plot
ivp.plotAllPumpProbe(path, full=full, 
                     autosave=autosave, autoclose=autoclose)

#%% LINEAR PREDICTION -------------------------------------------------------------

# Parameters
name = 'M_20190408_03.txt'
path = os.getcwd()
autoclose = True
autosave = True
round_Matlab_needed = True # Pyhon 3.6.2 needs it

# Load data
filename = os.path.join(path, name)
t, V, meanV, details = ivs.loadNicePumpProbe(filename)
t0 = ivp.interactiveTimeZero(filename, autoclose=autoclose)
t, V, meanV = iva.cropData(t0, t, V, meanV)
#t, V, meanV, details = iva.loadZeroPumpProbe(filename)
dt = details['dt']

# Use linear prediction
results, others = iva.linearPrediction(t, meanV, dt, 
                                      autoclose=autoclose,
                                      round_Matlab_needed=round_Matlab_needed)

# Plot linear prediction
ivp.linearPredictionPlot(filename, others, autosave=autosave)