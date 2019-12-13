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
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
path = os.path.join(home, r'Muestras\SEM\LIGO5bis\1')
series = 'LIGO5bis_1'

# Load data
rwidth = []
rheight = []
height = []
width = []
hangle = []
wangle = []
for file in os.listdir(path):
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

# Organize angle data...

# ...1st fix the horizontal angles measured upside down
new_hangle = []
for ha in hangle:
    new_ha = []
    for i in ha:
        difference = i - np.mean(ha)
        if abs(difference)>90:
            if abs(difference-180) < abs(difference+180):
                new_ha.append(i-180)
            else:
                new_ha.append(i+180)
        else:
            new_ha.append(i)
    new_hangle.append(new_ha)
    del new_ha, i
hangle = np.array(new_hangle).T
del new_hangle

# ...2nd fix the vertical angles measured upside down
new_wangle = []
for wa in wangle:
    new_wa = []
    for j in wa:
        difference = np.mean(wa) - j
        if abs(difference)>90:
            if abs(difference-180) < abs(difference+180):
                new_wa.append(j-180)
            else:
                new_wa.append(j+180)
        else:
            new_wa.append(j)
    new_wangle.append(new_wa)
    del new_wa, j
wangle = np.array(new_wangle).T
del new_wangle

# ...3rd rotate vertical angles to be horizontal ones
new_wangle = []
for ha, wa in zip(hangle.T, wangle.T):
    difference = np.mean(ha) - np.mean(wa)
    if abs(difference-90) < abs(difference+90):
        new_wangle.append(wa + 90)
    else:
        new_wangle.append(wa - 90)
wangle = np.array(new_wangle).T
del ha, wa, difference, new_wangle

# ...4th make all angles point between 0 and 135
angle = np.array([[*ha, *wa] for ha, wa in zip(hangle.T, wangle.T)]).T
new_angle = []
for a in angle.T:
    if np.mean(a) < 0:
        new_angle.append(a + np.ones(len(a))*180)
    elif np.mean(a) > 180:
        new_angle.append(a - np.ones(len(a))*180)
    else:
        new_angle.append(a)
angle = np.array(new_angle).T
del wangle, hangle, new_angle

# Get results
W = np.mean(width, axis=0)
dW = np.std(width, axis=0)
H = np.mean(height, axis=0)
dH = np.std(height, axis=0)
a = np.mean(angle, axis=0)
da = np.std(angle, axis=0)

# Apply correction due to method
H = H + dH
W = W + dW
A = H/W
dA = H*dW/W**2 + dH/W

# Organize results
results = np.array([W,dW,H,dH,A,dA,a,da]).T
heading = ["Ancho (nm)", "Error (nm)",
           "Longitud (nm)", "Error (nm)",
           "Relación de aspecto", "Error",
           "Ángulo (°)", "Error (°)"]

# Save data
ivs.saveTxt(
    os.path.join(path,'Resultados_SEM_{}.txt'.format(series)), 
    results,
    header=heading, footer=dict(rods=rods),
    overwrite=True
    )

# Round and gather results
items = []
for i in range(len(rods)):
    w = '\t'.join(ivu.errorValue(W[i], dW[i]))
    h = '\t'.join(ivu.errorValue(H[i], dH[i]))
    ra = '\t'.join(ivu.errorValue(A[i], dA[i], one_point_scale=True))
    an = '\t'.join(ivu.errorValue(a[i], da[i]))
    items.append('\t'.join([w, h, ra, an]))
del w, h, ra, an, W, H, A, a, dW, dH, dA, da

# Make OneNote table
heading = '\t'.join(heading)
items = ['\t'.join([n, r]) for n, r in zip(rods, items)]
items = '\n'.join(items)
heading = '\t'.join(['Rod', heading])
table = '\n'.join([heading, items])
ivu.copy(table)
del heading, items