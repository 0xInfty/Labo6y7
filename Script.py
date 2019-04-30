# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import os

#%% PARAMETERS --------------------------------------------------------------------

name = 'M_20190408_03.txt'
path = os.getcwd()
autoclose = True
autosave = True
round_Matlab_needed = True # Pyhon 3.6.2 needs it

#max_nrepetitions = 3 # BEWARE OF FALSE COLUMNS!

#%% LOAD AND CROP DATA ------------------------------------------------------------

# Load data
filename = os.path.join(path, name)
t, V, meanV, details = ivs.loadZeroPumpProbe(filename)
dt = details['dt']

#%% LINEAR PREDICTION -------------------------------------------------------------

# Use linear prediction
results, others = iva.linearPrediction(t, meanV, dt, 
                                      autoclose=autoclose,
                                      round_Matlab_needed=round_Matlab_needed)

#%% Plot linear prediction
ivp.linearPredictionPlot(filename, others, autosave=autosave)