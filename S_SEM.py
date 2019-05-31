# -*- coding: utf-8 -*-
"""
Created on Sat May 25 09:09:38 2019

@author: Vall
"""

import iv_utilities_module as ivu
import numpy as np
import os

path = r'C:\Users\quimica\Documents\Laboratorio\Profesores\Valeria Pais\Facu\OneDrive\Labo 6 y 7\SEM\M135\M135 Geometrías\5'
series = 'M135_5_1D'

names = []
height = []
width = []
for file in os.listdir(path):
#    print(file)
    if file.endswith("W.csv"):
        names.append(file.split('_W.csv')[0].split('_')[-1])
        width.append(np.loadtxt(os.path.join(path, file), 
                                delimiter=',', 
                                skiprows=1)[:,-1])
    elif file.endswith("H.csv"):
        height.append(np.loadtxt(os.path.join(path, file), 
                                 delimiter=',', 
                                 skiprows=1)[:,-1])
height = np.array(height).T
width = np.array(width).T
del file

W = np.mean(width, axis=0)
dW = np.std(width, axis=0)
H = np.mean(height, axis=0)
dH = np.std(height, axis=0)
A = W/H
dA = W*dH/H**2 + dW/H

results = np.array([W,dW,H,dH,A,dA]).T
heading = ''.join(["Ancho (nm)\tError (nm)\tAltura (nm)\t",
                   "Error (nm)\tRelación de aspecto\tError"])
np.savetxt(
    os.path.join(path,'Resultados_{}.txt'.format(series)), 
    results, 
    footer='\t'.join(names), 
    header=heading)

items = ['\t'.join([str(element) for element in row]) for row in results]
items = ['\t'.join([n, r]) for n, r in zip(names, items)]
items = '\n'.join(items)
heading = '\t'.join(['Nro', heading])
table = '\n'.join([heading, items])
ivu.copy(table)
