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
path = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7\Muestras\SEM\LIGO1\LIGO1 Geometrías\1'
series = 'LIGO1_1'

# Load data
rwidth = []
rheight = []
height = []
width = []
hangle = []
wangle = []
for file in os.listdir(path):
#    print(file)
    if file.endswith("W.csv"):
        rwidth.append(file.split('_W.csv')[0].split('_')[-1])
        width.append(np.loadtxt(os.path.join(path, file), 
                                delimiter=',', 
                                skiprows=1)[:,-1])
        wangle.append(np.loadtxt(os.path.join(path, file), 
                                 delimiter=',', 
                                 skiprows=1)[:,-2])
    elif file.endswith("H.csv"):
        rheight.append(file.split('_H.csv')[0].split('_')[-1])
        height.append(np.loadtxt(os.path.join(path, file), 
                                 delimiter=',', 
                                 skiprows=1)[:,-1])
        hangle.append(np.loadtxt(os.path.join(path, file), 
                                 delimiter=',', 
                                 skiprows=1)[:,-2])

# Organize length data
if rwidth!=rheight:
    raise ValueError("¡Falta algún dato!")
rods = rwidth
height = np.array(height).T
width = np.array(width).T
del file, rwidth, rheight

# Organize angle data
for ha, wa in hangle, wangle:
    difference = np.mean(ha) - np.mean(wa)
    standard = [90,-90]
    if abs(difference-standard[0]) < abs(difference-standard[1]):
        

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
