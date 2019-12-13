# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 09:48:07 2019

@author: Valeria
"""

from itertools import combinations as comb
import iv_save_module as ivs
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as st

#%% LINEAR PREDICTION REPRISE

# Filenames' Parameters
name = 'M_20191018_11'
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
series = 'Power' # 'CombinationsOf3' # 'EachOneByItsOwn'
path = os.path.join(home, 'Análisis', series + '_' + name)

# Experiments' Parameters
nexperiments = 4
groups_mode = 'own' # Combinations 'comb', each experiment by its own 'own'

# Analysis' Parameters
desired_frequency = 18 # GHz
max_frequency_deviation = 3 # GHz
nbins = 2

if False: print("""CombinationsOf3 / EachOneByItsOwn
                M_20191129...  01) 17.5 GHz\
                               02) 14.5 GHz\
                               03) 13 GHz""")

#%% MORE PARAMETERS

# Define other parameters
filenames = []
for file in os.listdir(os.path.join(path, 'Ajustes')):
    if 'Results' not in file and 'Figuras' not in file:
        filenames.append(os.path.join(path, 'Ajustes', file))

# Make groups of experiments
if groups_mode == 'comb':
    # Make all combinations of 3 elements
    experiments_groups = [list(c) for c in comb(list(range(nexperiments)), 3)]
elif groups_mode == 'own':
    # Take each experiment on its own
    experiments_groups = [[i] for i in range(nexperiments)]

#%% LOAD DATA

# Define variables and begin loop
all_groups = []
all_results = []
all_other_results = []
all_fit_params = []
#all_tables = []
for file in filenames:
    
    # Load data from a base fit made by hand
    results, header, footer = ivs.loadTxt(file)
    
    # Reorganize data
    other_results_keys = ['Nsingular_values', 'chi_squared']
    other_results = {k: footer[k] for k in other_results_keys}
    fit_params = dict(footer)
    for k in other_results_keys:
        fit_params.pop(k)
    del k, other_results_keys
    fit_params = ivu.InstancesDict(fit_params)
    del footer
    
#    # Generate fit tables
#    tables = iva.linearPredictionTables(fit_params, 
#                                        results, 
#                                        other_results))
    
    # Add data to external variables
    print(fit_params.use_experiments)
    all_groups.append(fit_params.use_experiments)
    all_results.append(results)
    all_other_results.append(other_results)
    all_fit_params.append(fit_params)
#    all_tables.append(tables)

del results, other_results, fit_params#, tables

# Make list of failed groups too
all_failed_groups = []
for g in experiments_groups:
    if g not in all_groups:
        all_failed_groups.append(g)
del g

#%% ANALYSIS

# Filter results and keep only desired frequencies
all_data = []
groups = []
results = []
other_results = []
fit_params = []
#tables = []
for i, g, fit in zip(range(len(all_groups)), all_groups, all_results):
        try:
            i = np.argmin(abs(fit[:,0] - desired_frequency*np.ones(fit.shape[0])))
            if abs(fit[i,0] - desired_frequency) > max_frequency_deviation:
                print("Group {} doesn't contain desired frequency".format(g))
                print(fit)
                all_failed_groups.append(g)
            else:
                all_data.append([*fit[i,:]])
                groups.append(g)
                results.append(all_results[i])
                other_results.append(all_other_results[i])
                fit_params.append(all_fit_params[i])
        except IndexError:
            if abs(fit[0] - desired_frequency) <= max_frequency_deviation:
                all_data.append(fit)
                groups.append(g)
                results.append(all_results[i])
                other_results.append(all_other_results[i])
                fit_params.append(all_fit_params[i])
#                tables.append(all_tables[i])
            else:
                print("Group {} doesn't contain desired frequency".format(g))
                print(fit)
                all_failed_groups.append(g)
del i, g, fit
all_data = np.array(all_data)
all_groups = groups
all_results = results
all_other_results = other_results
all_fit_params = fit_params
#all_tables = tables
del groups, results, other_results, fit_params#, tables

# Get interesting data
frequency = all_data[:,0]*1e9 # Hz
damping_time = all_data[:,1]*1e-12 # s
quality_factor = all_data[:,2]
chi_squared = [all_other_results[i]['chi_squared'] 
               for i in range(len(all_results))]
del all_other_results

# Save data
filenameMaker = lambda file : os.path.join(path, file)
new_header = list(header)
new_header.append('Chi cuadrado (V^2)')
ivs.saveTxt(filenameMaker('Results.txt'),
            np.array([*all_data.T, chi_squared]).T,
            header=new_header,
            footer={'experiments_groups': all_groups},
            overwrite=True)

# Print results
print('Serie ' + series + '\n')
print('Medición: ' + name)
print('Grupos de experimentos ajustados: {}/{}'.format(
        len(all_groups),
        len(experiments_groups)))
print(r"Frecuencia: {}".format(
        ivu.errorValueLatex(np.mean(frequency), 
                            np.std(frequency), 
                            symbol='±',
                            units="Hz")))
print(r"Tiempo de decaimiento: {}".format(
        ivu.errorValueLatex(np.mean(damping_time), 
                            np.std(damping_time), 
                            symbol='±',
                            units="s")))
print(r"Chi cuadrado: {}".format(
        ivu.errorValueLatex(np.mean(chi_squared),
                            np.std(chi_squared),
                            symbol='±',
                            units="V²")))

#%% HISTOGRAM WITH A CURVE

# Plot parameters
color = 'blue'
index = frequency.argsort()

# Make histogram
n, bins, patches = plt.hist(frequency*1e-9, nbins, density=True,
                            alpha=0.4, facecolor=color)
del patches

# Add curve over it
x = np.linspace(np.min(bins), np.max(bins), 50)
plt.plot(x, 
         st.norm.pdf(x, 
                     np.mean(frequency)*1e-9, 
                     np.std(frequency)*1e-9),
         color=color)
del index

# Format plot
plt.title('{} ({})'.format(series, name))
plt.xlabel("Frecuencia F (GHz)")
plt.ylabel(r"Densidad de probabilidad $\int f(F) dF = 1$")

# Save plot
ivs.saveFig(filenameMaker('Histogram.png'), overwrite=True)

if False: print(
""" ABOUT HISTOGRAM'S NORMALIZATION
- If either normed=1 or density=True, then it's True that...
  sum([n[i] * (bins[i+1]-bins[i]) for i in range(len(n))]) == 1
  ...meaning the area below the curve is equal to 1
- By default, then it's True that...
  sum(n) == len(n)
""")

if False: print(
""" ABOUT NORMAL DISTRIBUTION'S NORMALIZATION
By default, then the area below the curve is equal to 1.
""")

#%% BOXPLOT

fig = plt.figure()
grid = plt.GridSpec(1, 2, wspace=0.5, hspace=0)

ax = plt.subplot(grid[0,0])
ax.boxplot(frequency*1e-9, showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2})
for l in ax.get_xticklabels():
    l.set_visible(False)
del l
plt.ylabel("Frecuencia (GHz)")
ax.tick_params(axis='y', direction='in')

ax = plt.subplot(grid[0,1])
ax.boxplot(damping_time*1e12, showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2})
for l in ax.get_xticklabels():
    l.set_visible(False)
del l
plt.ylabel("Tiempo característico (ps)")
ax.tick_params(axis='y', direction='in')

# Save plot
ivs.saveFig(filenameMaker('Boxplots.png'), overwrite=True)

#%% PLOT

# Plot results for the different rods
fig, ax1 = plt.subplots()
plt.title('{} ({})'.format(series, name))

# Frequency plot, right axis
ax1.set_xlabel('Repetición')
ax1.set_ylabel('Frecuencia (GHz)', color='tab:red')
ax1.plot(frequency*1e-9, 'ro')
ax1.tick_params(axis='y', labelcolor='tab:red')

# Quality factor, left axis
ax2 = ax1.twinx()  # Second axes that shares the same x-axis
ax2.set_ylabel('Tiempo de decaimiento (ps)', color='tab:blue')
ax2.plot(damping_time*1e12, 'bx')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# Format graph
plt.grid(which='both', axis='x')
ax1.tick_params(length=2)
ax1.grid(axis='x', which='both')

# Save plot
ivs.saveFig(filenameMaker('FyTau.png'), overwrite=True)