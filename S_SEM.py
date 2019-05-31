# -*- coding: utf-8 -*-
"""
Created on Sat May 25 09:09:38 2019

@author: Vall
"""

import iv_utilities_module as ivu
import iv_save_module as ivs
import numpy as np
import os

# Parameters
path = r'C:\Users\quimica\Documents\Laboratorio\Profesores\Valeria Pais\Facu\OneDrive\Labo 6 y 7\SEM\M135\M135 Geometrías\7B\2D'
series = 'M135_7B_2D'

# Load data
rwidth = []
rheight = []
height = []
width = []
for file in os.listdir(path):
#    print(file)
    if file.endswith("W.csv"):
        rwidth.append(file.split('_W.csv')[0].split('_')[-1])
        width.append(np.loadtxt(os.path.join(path, file), 
                                delimiter=',', 
                                skiprows=1)[:,-1])
    elif file.endswith("H.csv"):
        rheight.append(file.split('_H.csv')[0].split('_')[-1])
        height.append(np.loadtxt(os.path.join(path, file), 
                                 delimiter=',', 
                                 skiprows=1)[:,-1])

# Organize data
if rwidth!=rheight:
    raise ValueError("¡Falta algún dato!")
rods = rwidth
height = np.array(height).T
width = np.array(width).T
del file, rwidth, rheight

# Get results
W = np.mean(width, axis=0)
dW = np.std(width, axis=0)
H = np.mean(height, axis=0)
dH = np.std(height, axis=0)
A = W/H
dA = W*dH/H**2 + dW/H

# Organize results
results = np.array([W,dW,H,dH,A,dA]).T
heading = ["Ancho (nm)", "Error (nm)",
           "Altura (nm)", "Error (nm)",
           "Relación de aspecto", "Error"]

# Save data
ivs.saveTxt(
    os.path.join(path,'Resultados_{}.txt'.format(series)), 
    results,
    header=heading,
    footer=dict(rods=rods)
    )

# Make OneNote table
heading = ''.join(heading)
items = ['\t'.join([str(element) for element in row]) for row in results]
items = ['\t'.join([n, r]) for n, r in zip(rods, items)]
items = '\n'.join(items)
heading = '\t'.join(['Rod', heading])
table = '\n'.join([heading, items])
ivu.copy(table)
