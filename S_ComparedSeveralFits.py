# -*- coding: utf-8 -*-
"""
Created on Mon Dec 09 23:28:53 2019

@author: Vall
"""

from itertools import combinations as comb
import iv_save_module as ivs
import iv_utilities_module as ivu
import numpy as np
import os

#%% LINEAR PREDICTION REPRISE

# Filenames' Parameters
names = ['M_20191129_01', 'M_20191129_02', 'M_20191129_03']
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
nexperiments = [10, 10, 10]

# Groups' Parameters
groups_modes = ['own', 'comb'] # Combinations 'comb', each experiment by its own 'own'
series = ['EachOneByItsOwn', 'CombinationsOf3']

# Analysis' Parameters
desired_frequency = [17.5, 14.5, 13] # GHz
max_frequency_deviation = [2, 2, 2] # GHz

# Save parameters
paths = [[os.path.join(home, 'Análisis', s + '_' + n) for n in names] for s in series]
overwrite = False
autosave = True

#%% MORE PARAMETERS

# Make groups of experiments
experiments_groups = []
for i, g in enumerate(groups_modes):
    if g == 'comb':
        # Make all combinations of 3 elements
        experiments_groups.append(
                [list(c) for c in comb(list(range(nexperiments[i])), 3)])
    elif g == 'own':
        # Take each experiment on its own
        experiments_groups.append([[i] for i in range(nexperiments[i])])
del i, g

#%% LOAD DATA

# Define variables and begin loop
all_groups = []
all_data = []
for i, pi in enumerate(paths):
    these_groups = []
    this_data = []
    for j, pij in enumerate(pi):
        # Load data
        results, header, footer = ivs.loadTxt(os.path.join(pij, 'Results.txt'))        
        # Export data outside loop
        these_groups.append(footer['experiments_groups'])
        this_data.append(results)
    all_groups.append(these_groups)
    all_data.append(this_data)
del i, pi, footer, results, these_groups, this_data

'''
# Make list of failed groups too
all_failed_groups = []
for groups, egroups in zip(all_groups, experiments_groups):
    failed_groups = []
    for g in groups:
        for eg in egroups:
            if g not in groups:
                print(str(g) + '\n')
                print(str(groups) + '\n\n')
                failed_groups.append(g)
        all_failed_groups.append(failed_groups)
del groups, failed_groups, g, egroups

#%% ANALYSIS

# Get interesting data
frequency = [all_data[i][:,0]*1e9 # Hz
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
                            units="Pa")))
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

#%% MAKE HISTOGRAM WITH A CURVE

# Plot parameters
nbins = 20
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

#%% MAKE BOXPLOT

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
'''