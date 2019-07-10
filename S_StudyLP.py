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

#%%

# Parameters
home = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7'
name = 'M_20190605_11'

# Save parameters
autosave = True
overwrite = True

# Plot parameters
plot_params = dict(
        plot = False,
        interactive = False,
        autoclose = True,
        )
plot_params = ivu.InstancesDict(plot_params)

# Fit parameters
fit_params = dict(
        use_full_mean = True,
        use_experiments = [1], # First is 0, not 1!
        send_tail_to_zero = False,
        tail_method = 'mean', # Could also be 'min' or 'max' or any numpy function
        use_fraction = .2,
        choose_t0 = True,
        choose_tf = False,
        max_svalues = 8,
        )
fit_params = ivu.InstancesDict(fit_params)

# Create full filename
filename = ivs.filenameToMeasureFilename(name, home)

# Load data
t, V, details = ivs.loadNicePumpProbe(filename)

# Choose data to fit
if fit_params.use_full_mean:
    data = np.mean(V, axis=1)
else:
    data = np.mean(V[:, fit_params.use_experiments], axis=1)

# Choose time interval to fit
t0 = -8.122954823529412 # This is an initial time we think that optimizes it
i = np.argmin(np.abs(t-t0)) # We'll take this index as main initial time
Ni = 40 # We'll try this many index to the right and left from the main index

#%%

# Now iterate, fitting on different initial times
results = []
other_results = []
fit_terms = []
jgood = [] # Here we'll collect the index that allow fitting
for j in range(max(i-Ni,0), i+Ni):

    # Choose initial time t0
    t0j = t[j]
#    print(t0j, j)
    
    # Crop data accorddingly
    tj, dataj = iva.cropData(t0j, t, data)
    fit_params.time_range = (t0j, t[-1])
    fit_params.voltage_zero = 0
    
    # Use linear prediction, if allowed
    try:
        res, other, plot = iva.linearPrediction(tj, 
                                                dataj, 
                                                details['dt'], 
                                                svalues=4,
                                                printing=False)
        jgood.append(j)
        results.append(res)
        other_results.append(other)
        fit_terms.append(plot.fit)
    except:
        results.append(None)
        other_results.append(None)
        fit_terms.append(None)
        
del t0j, tj, dataj, res, other, plot

# Now select only the fits that satisfy us
jreallygood = []
jrare = [] # The ones that hold only one oscillant term
frequencies = []
quality = []
chi = []
meanqdiff = []
stdqdiff = []
for j in jgood:
    res = results[j]
    other = other_results[j]
    if res.shape[0]!=1:
        imax = np.argmin(np.abs(res[:,0] - 9 * np.ones(len(res[:,0]))))
        if res[imax,0] != 0:
            frequencies.append(res[imax,0])
            quality.append(res[imax,2])
            chi.append(other['chi_squared'])
            jreallygood.append(j)
            term = fit_terms[j]
            meanqdiff.append( np.mean( (term[:,1]-term[:,imax+2])**2 ) )
            stdqdiff.append( np.std( (term[:,1]-term[:,imax+2])**2 ))
    else:
        if res[0,0] != 0:
            frequencies.append(res[0,0])
            quality.append(res[0,2])
            chi.append(other['chi_squared'])
            jreallygood.append(j)
            jrare.append(j)
            term = fit_terms[j]
            meanqdiff.append( np.mean( (term[:,1]-term[:,3])**2 ) )
            stdqdiff.append( np.std( (term[:,1]-term[:,imax+2])**2 ))
del res, other

#%%

plt.figure()
plt.plot(jreallygood, frequencies, 'x')
plt.plot(i, frequencies[i], 'xr')
plt.ylabel('Frecuencia (GHz)')

plt.figure()
plt.plot(jreallygood, quality, 'o')
plt.plot(i, quality[i], 'or')
plt.ylabel('Factor de calidad')

plt.figure()
plt.plot(jreallygood, chi, '.')
plt.plot(i, chi[i], 'xr')
plt.ylabel('Chi cuadrado')

plt.figure()
plt.plot(jreallygood, meanqdiff, 'x')
plt.plot(i, meanqdiff[i], 'xr')
plt.ylabel('Diferencia cuadrática media')

plt.figure()
plt.plot(jreallygood, stdqdiff, 'x')
plt.plot(i, stdqdiff[i], 'xr')
plt.ylabel('Desviación estándar de la diferencia cuadrática')