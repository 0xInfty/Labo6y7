# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 10:57:00 2019

@author: Valeria
"""

import iv_analysis_module as iva
import iv_save_module as ivs
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import os

#%% PARAMETERS

home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
folders = [r'M135\7B\1', r'M135\5\1']
series = ['M135_7B_1', 'M135_5_1']
full_series = ['SiO$_2$ (Grilla 7B)', r'SiO$_2$ (Grilla 5)']

#home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
#folders = [r'LIGO1\1', r'LIGO5bis\1']
#series = ['LIGO1_1', 'LIGO5bis_1']
#full_series = ['L1 (FS)', r'L5 (Ta$_2$O$_5$)']

filter_ta2o5_outliers = True # largest L

figsFilename = lambda n : os.path.join(home, r'Análisis\LoadSEM', n+'.png')
figs_extension = '.png'
symbols=['% ', '$\pm$'] #['', '±']
make_boxplot_of = [[0], [1], [0,1]]
overwrite = True

#%% LOAD DATA

# Organize paths
paths = [os.path.join(home, r'Muestras\SEM', f) for f in folders]
filenames = [os.path.join(home, r'Muestras\SEM', 
                          f, 'Resultados_SEM_{}.txt'.format(s))
             for f, s in zip(folders, series)]

data = []
rods = []
for f in filenames:
    d, header, ft = ivs.loadTxt(f)
    data.append(d)
    rods.append(ft['rods'])
del d, ft

if filter_ta2o5_outliers:
    index = np.argsort(data[1][:,2])[:-1]
    rods[1] = [rods[1][i] for i in index]
    data[1] = data[1][index,:]

#%% VALUES

variables = ['Longitud L', 'Diámetro d', 'Relación de aspecto', 'Ángulo']
variables_units = ['nm', 'nm', '', 'º']
variables_data = lambda i : [
        iva.getValueError(data[i][:,2], data[i][:,3]),
        iva.getValueError(data[i][:,0], data[i][:,1]),
        iva.getValueError(data[i][:,4], data[i][:,5]),
        iva.getValueError(data[i][:,6], data[i][:,7])]

values_string = '{}Serie LoadSEM\n\n'.format(symbols[0])
for i, fs in enumerate(full_series):
    values_string += "{}Resultados de {}\n".format(symbols[0], fs)
    values_string += "{}Cantidad de NPs: {:.0f}\n".format(symbols[0],
                                                          len(data[i][:,0]))
    for v, vd, vu in zip(variables, variables_data(i), variables_units):
        values_string += '{}{} = {}\n'.format(symbols[0], 
                                              v, 
                                              ivu.errorValueLatex(
                                                      *vd, 
                                                      units=vu, 
                                                      symbol=symbols[1]))
    values_string += '\n'
print(values_string)
ivu.copy(values_string)

#%% *) BOXPLOTS L, d

for mbo in make_boxplot_of:
    population = [data[s].shape[0] for s in range(len(full_series))]
    labels = [('{}\n'+r'$\hookrightarrow${} NPS').format(full_series[s], population[s]) 
              for s in mbo]  
    
    # Data inside each series to plot
    choose_index_from_header = [2, 0]#, 4]
    boxplot_data = [[data[j][:,i] for j in mbo] 
                    for i in choose_index_from_header]
    ax_labels = [r"Longitud $L$ (nm)", r"Diámetro $d$ (nm)"]#, r"Ángulo $\Phi$"]
    
    # Format
    base_height = .1
    base_width = .6
    label_left_space = .08
    label_right_space = .08
    if len(mbo)==1:
        alpha = [.25, .5, .5]
    elif len(mbo)==2:
        alpha = [0, .15, 0]
    elif len(mbo)==3:
        alpha = [0, .05, 0]
    
    # Begin Figure
    fig = plt.figure()
    grid = plt.GridSpec(len(choose_index_from_header), 1, hspace=0.1)
    ax = [plt.subplot(g) for g in grid]
    
    index = 0
    for a, dat, lab in zip(ax, boxplot_data, ax_labels):
        
        # Boxplot
        bplot = a.boxplot(
            dat, 
            showmeans=True, meanline=True, 
            meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
            medianprops={'color':'r', 'linewidth':2},
            flierprops={'markersize':7},
            patch_artist=True,
            widths=base_width,
            labels=labels,
            vert=False)
        for p in bplot['boxes']:
            p.set_facecolor('w') # paint white boxes
        del p, bplot    
    
        # Labels' format
        a.xaxis.set_label_text(lab, va='center')
    #    ax.tick_params(axis='x', direction='in')
    
        # Grid's format
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.grid(which='major', axis='x')
        a.grid(which='minor', axis='x', linestyle=':')
        a.grid(which='major', axis='y')
        a.yaxis.tick_right()
        a.yaxis.set_label_position('right')
        
        # Axes size
        box = a.get_position()
    #    if len(mbo)!=1:
        w = box.width
        box.x0 = box.x0 - w * label_left_space
        box.x1 = box.x1 - w * label_right_space
        box.y1 = box.y0 + base_height * len(mbo)
        box.y0 = box.y0 + alpha[index]
        box.y1 = box.y1 + alpha[index]
        a.set_position(box)
        
    #    # Add population per box
    #    population = [len(d) for d in dat]
    #    positions = list(a.get_yticks())
    #    for n, pos in zip(population, positions):
    #        a.text(-.1, pos, '{:.0f} NPs'.format(n))
        
        index += 1
    
    del index, box, w, a
    
    ax[0].xaxis.tick_top()
    ax[0].xaxis.set_label_position('top')
    
    #fig.text(.9, .9, '{')
    
    ivs.saveFig(figsFilename('Boxplots{}'.format(mbo)), 
                overwrite=overwrite)
    
    del mbo, population, labels
    del base_height, base_width, label_right_space, label_left_space
    del choose_index_from_header, boxplot_data, ax_labels
    del ax, grid
    
