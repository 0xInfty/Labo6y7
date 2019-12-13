# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 09:48:07 2019

@author: Vall
"""

import iv_analysis_module as iva
import iv_save_module as ivs
import iv_utilities_module as ivu

#%% PARAMETERS -------------------------------------------------------------------

# Parameters
fit_name = 'M_20191129_01'
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'

#%% LOAD FIT DATA

# Create fit filename
fit_filename = ivs.filenameToFitsFilename(fit_name, home=home)

# Load data
results, header, footer = ivs.loadTxt(fit_filename)

# Reorganize data
other_results_keys = ['Nsingular_values', 'chi_squared']
other_results = {k: footer[k] for k in other_results_keys}
fit_params = dict(footer)
for k in other_results_keys:
    fit_params.pop(k)
del k, other_results_keys
fit_params = ivu.InstancesDict(fit_params)
del footer

# Generate tables
tables = iva.linearPredictionTables(fit_params, results, other_results)
ivu.copy(tables[0])